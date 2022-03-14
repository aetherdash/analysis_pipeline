import os, sys
import random
import pandas as pd
import numpy as np
import pickle
import json
import itertools
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

# Classification models
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.neural_network import MLPClassifier

# Regression models
from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.neural_network import MLPRegressor

from analytics_utils.database_access.table_properties import *
from analytics_utils.database_access.sql_utils import get_metric_filter_str
from analytics_utils.database_access.db_interface import DatabaseInterface
from analytics_utils.database_access.s3_interface import download_from_s3, upload_to_s3, s3_imgupload, s3_df2csv
from analytics_utils.maldi_processing.ttest import get_peak_vals, get_ttests
from analytics_utils.analysis_tools.dataset_creation_utils import get_source_plate_grp, update_source_grp, combine_rows_by_unique_enzyme, combine_maldi_lcms_data, get_ensemble_cv_results, get_ensemble_predictions
from analytics_utils.analysis_tools.analysis_utils import get_ctrl_vals, get_lcms_qc_derivedmetrics_checks, calculate_conversion_enantioselectivity_C18, calculate_conversion_enantioselectivity_chiral
from analytics_utils.utils.slack_interface import post_slack_message, SlackMessageConstants
from analytics_utils.visualization_tools.visualization_utils import generate_plots, df2array_dict
from analytics_utils.lims_tools.lims_utils import lims_post_matplotlib_img

        
class DatasetCreation:
    
    def __init__(self,
                 s3_bucket='ml-analytics-file-store',
                 s3_subfolder='ProjectX_UnitX/',
                 model_fname='model',
                 fname_suffix='',
                 maldi_detections_1='maldi_detections_1',
                 maldi_detections_2='maldi_detections_2',
                 lcms_detections='lcms_detections',
                 lcms_detections_chiral='lcms_detections_chiral',
                 maldi_classification_features=maldi_classification_feature_names,
                 lcms_classification_features=lcms_classification_feature_names,
                 model_type_list=['maldi1sample_classification', 'maldi2sample_classification', 'maldi2lcms_regression'],
                 combine_maldi_lcms_by = ['source_address', 'source_plate_grp'],
                 posctrl_enzyme_barcode = None,
                 maldi_label_condition = {},
                 hit_filter_dict = {},
                 plot_enantiomer_lcms_metrics = None,
                 col_suffix_db = {'':'_(r)', 'CMP60354':'_(r)', 'CMP60403':'_(-)', 'CMP60404':'_(+)'},
                 neg_ctrltype='EV',
                 extra_peaks = 0,
                 normalization_peaks = [172,379],
                 local_dir='../DATASETS/',
                 update_plate_analytics_table=True,
                 update_enzyme_analytics_table=True,
                 track_dilution_plates=False,
                 save_plots_to_lims=True,
                 get_2sample_from_1sample_dataset=True
                ):

        self.s3_bucket = s3_bucket
        self.s3_subfolder = s3_subfolder
        self.model_fname = model_fname
        self.fname_suffix = fname_suffix
        self.model_path = f'{self.s3_subfolder}{self.model_fname}'
        self.local_dir = local_dir
        self.maldi_detections_1 = maldi_detections_1
        self.maldi_detections_2 = maldi_detections_2
        self.lcms_detections = lcms_detections
        self.lcms_detections_chiral = lcms_detections_chiral
        self.maldi_classification_features = maldi_classification_features
        self.lcms_classification_features = lcms_classification_features
        self.lcms_output_feature = 'prod_conc_lcms_actual'
        self.neg_ctrltype = neg_ctrltype
        self.extra_peaks = extra_peaks
        self.normalization_peaks = normalization_peaks
        self.update_plate_analytics_table = update_plate_analytics_table
        self.update_enzyme_analytics_table = update_enzyme_analytics_table
        self.track_dilution_plates = track_dilution_plates
        self.model_type_list = model_type_list
        self.message_suffix = ''
        self.cols_to_get = {'maldi_detections_1': maldi_detections_columns_to_get,
                            'maldi_detections_2': maldi_detections_columns_to_get,
                            'lcms_detections': lcms_detections_columns_to_get,
                            'lcms_detections_chiral': lcms_detections_chiral_columns_to_get}
        self.combine_maldi_lcms_by = combine_maldi_lcms_by
        self.posctrl_enzyme_barcode = posctrl_enzyme_barcode
        self.maldi_label_condition = maldi_label_condition
        self.hit_filter_dict = hit_filter_dict
        self.plot_enantiomer_lcms_metrics = plot_enantiomer_lcms_metrics
        self.col_suffix_db = col_suffix_db
        self.save_plots_to_lims = save_plots_to_lims
        self.get_2sample_from_1sample_dataset = get_2sample_from_1sample_dataset
        self.ctrl_to_color_mapping = {'exp':'r', 'pos':'g', 'EV':'k'}
        self.cols_to_consolidate = {
            'C18_chiral': ['dev_or_prod', 'exp_workflow_barcode', 'exp_workflow_name', 'proj_barcode', 'proj_name', 'ctrl_type', 'exp_condition', 'enzyme_barcode', 'sequence', 'hamming', 'mutations', 'reference_enzyme', 'enzyme_concentration', 'enzyme_unit', 'enzyme_class', 'sequence_qc'] + library_metadata + plate_time_metadata + pelletqc_metadata,
            'maldi_lcms': ['dev_or_prod', 'exp_workflow_barcode', 'exp_workflow_name', 'proj_barcode', 'proj_name', 'ctrl_type', 'exp_condition', 'enzyme_barcode', 'sequence', 'hamming', 'mutations', 'reference_enzyme', 'enzyme_concentration', 'enzyme_unit', 'enzyme_class', 'sequence_qc'] + library_metadata + plate_time_metadata + pelletqc_metadata
        }
    
    
    def load_data(self, data_type, plate=None, run=None, exp_workflow=None, csv_fname=None):
        """
        Load maldi or lcms dataset from Postgres or from local csv file.
        """
        if csv_fname is None:
            try:
                # load from Postgres table
                filter_dict = {}
                if plate is not None: filter_dict.update({'plate':plate})
                if run is not None: filter_dict.update({'run':run})
                if exp_workflow is not None: filter_dict.update({'exp_workflow_barcode': exp_workflow})
                if len(filter_dict) > 0:
                    df = DatabaseInterface(data_type).lookup_previous_record(filter_dict, self.cols_to_get[data_type])
            except:
                print(f'Unable to load {data_type} data.')
        else:
            try:
                # load from CSV file saved locally
                fname = f"{self.local_dir}{csv_fname}"
                print(fname)
                df = pd.read_csv(fname, index_col=0, thousands=',')
                print(f'Loaded {data_type} data.')
            except:
                print(f'Unable to load {data_type} data.')
        
        # remove any pos or exp ctrl wells in which enzyme_barcode == nan
        df = df[df.ctrl_type.isin([self.neg_ctrltype,'pos','exp'])]
        df.loc[df.ctrl_type==self.neg_ctrltype, 'enzyme_barcode'] = 'DNA10001'
        if self.posctrl_enzyme_barcode is not None:
            df.loc[df.ctrl_type=='pos', 'enzyme_barcode'] = self.posctrl_enzyme_barcode
        df = df[(df.ctrl_type==self.neg_ctrltype) | ((df.ctrl_type=='pos') & (~df.enzyme_barcode.isnull())) | (df.ctrl_type=='exp')]
            
        return df
    
    
    def get_maldi_1sample_data(self, plate, run, exp_workflow, maldi_1_csv,
                               get_source_plate_groupings=True,
                               source_plate_groupings={},
                               plate_to_grp_idx_mapping={}):
        """
        Get MALDI 1-sample dataset
        """
        df = self.load_data(self.maldi_detections_1, plate, run, exp_workflow, maldi_1_csv)
        print('Obtained MALDI 1-sample dataset.')
        self.plate_list = list(set(df.plate))
        source_plate_list = list(set(df.source_plate))
        self.substrate_list = {plate: list(set(df.loc[df.plate==plate, 'substrate_barcode'])) for plate in self.plate_list}
        if get_source_plate_groupings:
            self.source_plate_groupings, self.plate_to_grp_idx_mapping = get_source_plate_grp(df, plate_str_list=['rxn', 'main'],
                                                                                              source_plate_list=source_plate_list,
                                                                                              source_plate_groupings=source_plate_groupings,
                                                                                              plate_to_grp_idx_mapping=plate_to_grp_idx_mapping)
        
        df = update_source_grp(df, self.source_plate_groupings)
        return df


    def get_maldi_2sample_data(self, plate, run, exp_workflow, maldi_2_csv,
                               df_maldi_1=None, simulate_2sample_dataset=False,
                               num_ref=16, num_rxn=4, num_unique_samples=48,
                               get_source_plate_groupings=True,
                               source_plate_groupings={},
                               plate_to_grp_idx_mapping={}):
        """
        Get 2-sample maldi data
        """
        if simulate_2sample_dataset:
            df = self.simulate_2sample_from_1sample_dataset(df_maldi_1, num_ref, num_rxn, num_unique_samples)
        elif self.get_2sample_from_1sample_dataset:
            self.plate_list = list(set(df_maldi_1.plate))
            df = pd.DataFrame(columns=list(df_maldi_1.columns))
            for plate in self.plate_list:
                print(plate, end=' ')
                maldi1_plate = df_maldi_1[df_maldi_1.plate==plate]
                maldi2_plate = self.simulate_2sample_from_1sample_dataset(maldi1_plate, group_by='source_address')
                df = df.append(maldi2_plate, ignore_index=True)
        else:
            df = self.load_data(self.maldi_detections_2, plate, run, exp_workflow, maldi_2_csv)
        print('Obtained MALDI 2-sample dataset.')
        self.plate_list = list(set(df.plate))
        self.substrate_list = {plate: list(set(df.loc[df.plate==plate, 'substrate_barcode'])) for plate in self.plate_list}
        source_plate_list = list(set(df.source_plate))
        if get_source_plate_groupings:
            self.source_plate_groupings, self.plate_to_grp_idx_mapping = get_source_plate_grp(df, plate_str_list=['rxn', 'main'],
                                                                                              source_plate_list=source_plate_list,
                                                                                              source_plate_groupings=source_plate_groupings,
                                                                                              plate_to_grp_idx_mapping=plate_to_grp_idx_mapping)
        df = update_source_grp(df, self.source_plate_groupings)
        return df
    

    def get_lcms_C18_data(self, lcms_run, exp_workflow, lcms_csv, get_source_plate_groupings=True, source_plate_groupings={}, plate_to_grp_idx_mapping={}):
        df = self.load_data(self.lcms_detections, None, lcms_run, exp_workflow, lcms_csv)
        df = df[~df.plate.isnull()]
        source_plate_list = list(set(df.source_plate))
        if get_source_plate_groupings:
            self.source_plate_groupings, self.plate_to_grp_idx_mapping = get_source_plate_grp(df, plate_str_list=['rxn', 'main'],
                                                                                              source_plate_list=source_plate_list,
                                                                                              source_plate_groupings=source_plate_groupings,
                                                                                              plate_to_grp_idx_mapping=plate_to_grp_idx_mapping)
        df = update_source_grp(df, self.source_plate_groupings)
        print('Obtained LCMS C18 dataset.')
        return df
    
    def get_lcms_chiral_data(self, lcms_run, exp_workflow, lcms_csv, get_source_plate_groupings=True, source_plate_groupings={}, plate_to_grp_idx_mapping={}):
        df = self.load_data(self.lcms_detections_chiral, None, lcms_run, exp_workflow, lcms_csv)
        df = df[~df.plate.isnull()]
        source_plate_list = list(set(df.source_plate))
        if get_source_plate_groupings:
            self.source_plate_groupings, self.plate_to_grp_idx_mapping = get_source_plate_grp(df, plate_str_list=['rxn', 'main'],
                                                                                              source_plate_list=source_plate_list,
                                                                                              source_plate_groupings=source_plate_groupings,
                                                                                              plate_to_grp_idx_mapping=plate_to_grp_idx_mapping)
        df = update_source_grp(df, self.source_plate_groupings)
        print('Obtained LCMS chiral dataset.')
        return df
    
    
    def sample_replicates_dataset(self, ref_data, rxn_data, num_ref, num_rxn, num_unique_samples, replicates_data={'ref':[], 'rxn':[]}):
        """
        From 1-sample detections data, sample 2-sample replicate data with a given number of ref and rxn replicates.
        """
        from utils.additional_utils import split_list
        # split rxn data
        num_split_rounds_rxn = int(np.ceil(num_unique_samples / (len(rxn_data)/num_rxn)))
        for r in range(num_split_rounds_rxn):
            random.seed(r)
            rxn_idxs = list(range(len(rxn_data)))
            random.shuffle(rxn_idxs)
            rxn_idxs_split = split_list(rxn_idxs, chunk_size=num_rxn)
            rxn_data_split = [rxn_data.loc[rxn_idxs] for rxn_idxs in rxn_idxs_split]
            replicates_data['rxn'] += rxn_data_split

        # split ref_data
        num_unique_samples = len(replicates_data['rxn'])
        num_split_rounds_ref = int(np.ceil(num_unique_samples / (len(ref_data)/num_ref)))
        for r in range(num_split_rounds_ref):
            ref_idxs = list(range(len(ref_data)))
            random.shuffle(ref_idxs)
            ref_idxs_split = split_list(ref_idxs, chunk_size=num_ref)
            ref_data_split = [ref_data.loc[ref_idxs] for ref_idxs in ref_idxs_split]
            replicates_data['ref'] += ref_data_split

        return replicates_data
    
    def group_replicates_dataset(self, ref_data, rxn_data, group_by='source_address', replicates_data={'ref':[], 'rxn':[]}):
        """
        Group 1-sample detections data into sets of replicates, based on given column
        """
        # group exp and pos control data
        unique_groups_rxn = list(set(rxn_data[group_by]))
        for grp in unique_groups_rxn:
            rxn_data_split = rxn_data[rxn_data[group_by]==grp]
            replicates_data['rxn'].append(rxn_data_split)
        
        # add reference control groups to 'rxn' replicates
        unique_groups_ref = list(set(ref_data[group_by]))
        for grp in unique_groups_ref:
            ref_data_split = ref_data[ref_data[group_by]==grp]
            replicates_data['rxn'].append(ref_data_split)
                
        # get reference replicates
        replicates_data['ref'] += [ref_data]*len(replicates_data['rxn'])
        return replicates_data
    
            
    def simulate_2sample_from_1sample_dataset(self, df, num_ref=None, num_rxn=None, num_unique_samples=None, neg_ctrltype='EV', group_by='source_address'):
        """
        Get sampled replicates data, and calculate mean and T-stats on replicates
        """
        # parse dataframe by substrate barcode
        data_sorted = {}
        substrate_barcode_list = list(set(df.substrate_barcode))
        for substrate_barcode in substrate_barcode_list:
            # filter by substrate
            df_subset = df[df.substrate_barcode==substrate_barcode]
            # get 'ref' data
            ref_data = df_subset[df_subset.ctrl_type==neg_ctrltype]
            ref_data.reset_index(drop=True, inplace=True)
            # get 'rxn' data
            rxn_data = df_subset[df_subset.ctrl_type!=neg_ctrltype]
            rxn_data.reset_index(drop=True, inplace=True)
            data_sorted.update({substrate_barcode: {'ref': ref_data, 'rxn': rxn_data}})
        
        # get replicates data
        replicates_data = {'ref':[], 'rxn':[]}
        # randomly sample from all 1-sample data with sample sizes as defined in input
        if num_ref is not None and num_rxn is not None and num_unique_samples is not None:
            for substrate_barcode, data_dict_substrate in data_sorted.items():
                ref_data = data_dict_substrate['ref']
                rxn_data = data_dict_substrate['rxn']
                replicates_data = self.sample_replicates_dataset(ref_data, rxn_data, num_ref, num_rxn, num_unique_samples,
                                                                 replicates_data=replicates_data)
        # group replicates data according to common value in designated column
        else:
            for substrate_barcode, data_dict_substrate in data_sorted.items():
                ref_data = data_dict_substrate['ref']
                rxn_data = data_dict_substrate['rxn']
                replicates_data = self.group_replicates_dataset(ref_data, rxn_data, group_by='source_address',
                                                                 replicates_data=replicates_data)
            
        
        # get mean and T-statistics
        # initialize dataframe
        simulated_maldi_data_2 = pd.DataFrame(columns=list(replicates_data['ref'][0].columns))

        for i, (ref_data, rxn_data) in enumerate(zip(replicates_data['ref'], replicates_data['rxn'])):
            d = dict(rxn_data.iloc[0])
            for colname in ['address', 'id', 'spectra_ids', 'source_address']:
                vals = sorted(list(set(list(rxn_data[colname]))))
                d[colname] = ', '.join([str(val) for val in vals])
            d["spectra_ids_ref"] = ', '.join(list(ref_data["spectra_ids_ref"]))
            if rxn_data.iloc[0]['true_score_binary'] is None or np.isnan(rxn_data.iloc[0]['true_score_binary']):
                d['true_score_binary'] = np.nan
            else:
                d['true_score_binary'] = (float(rxn_data['true_score_binary'].mean()) > 0.5)*1
                
            # GET AVG PEAK VALS
            d = get_peak_vals(d, ref_data, rxn_data, extra_peaks=self.extra_peaks, normalization_peaks=self.normalization_peaks)
            
            # PERFORM T-TESTS
            d = get_ttests(d, ref_data, rxn_data, extra_peaks=self.extra_peaks, normalization_peaks=self.normalization_peaks)

            # APPEND TO DATAFRAME
            simulated_maldi_data_2 = simulated_maldi_data_2.append(pd.DataFrame(d, index=[i]))
        
        return simulated_maldi_data_2
    
    
    def label_maldi2_with_maldi1(self, df_2, df_1_labeled, label_col, spectra_ids_col, binarize=True):
        df_2_labeled = df_2.copy()
        df_2_labeled[label_col] = 0
        for i in list(df_2_labeled.index):
            spectra_ids = str(df_2_labeled.loc[i, spectra_ids_col]).split(', ')
            df_1_replicates = df_1_labeled[df_1_labeled[spectra_ids_col].isin(spectra_ids)]
            true_score = df_1_replicates[label_col].mean()
            if binarize:
                true_score = (true_score>=0.5)*1
            df_2_labeled.loc[i, label_col] = true_score
        return df_2_labeled
    
    
    def get_binary_labels_from_maldi_data(self, df_1, df_2, split_by_list=[''], neg_ctrltype='EV',
                                          binary_label_col='true_score_binary',
                                          maldi_label_condition={'pk_prod_stdz_379':{'pk_thres':0.045, 'rxn_thres_factor':1}}):
        
        for col_suffix in split_by_list:
            # get col names
            label_col = f'{binary_label_col}_{col_suffix}' if col_suffix!='' else binary_label_col
            
            # label 1-sample maldi data by maldi peaks
            df_1_labeled = df_1.copy()
            df_1_labeled[label_col] = 0
            for metric, conditions in maldi_label_condition.items():
                metric_col = f'{metric}_{col_suffix}' if col_suffix!='' else metric
                pk_thres = conditions['pk_thres']
                rxn_thres_factor = conditions['rxn_thres_factor']
                # label neg ctrl data above pk_thres as '1'
                df_1_labeled.loc[((df_1_labeled['ctrl_type']==neg_ctrltype) & (df_1_labeled[metric_col]>pk_thres)), label_col] = 1
                # label pos and exp ctrl data above pk_thres * rxn_thres_factor as '0'
                df_1_labeled.loc[(df_1_labeled['ctrl_type'].isin(['exp','pos']) & (df_1_labeled[metric_col]>pk_thres*rxn_thres_factor)), label_col] = 1

            for ctrl_type in [self.neg_ctrltype, 'pos', 'exp']:
                df_1_labeled_ctrltype = df_1_labeled[df_1_labeled.ctrl_type==ctrl_type]
                print(f"[{ctrl_type}] Hit rate: {round(df_1_labeled_ctrltype[label_col].mean(),4)} (n={len(df_1_labeled_ctrltype)})")

            # label 2-sample maldi data by labeled 1-sample maldi data
            if df_2 is not None:
                spectra_ids_col = f'spectra_ids_{col_suffix}' if col_suffix!='' else 'spectra_ids'
                df_2_labeled = self.label_maldi2_with_maldi1(df_2, df_1_labeled, label_col, spectra_ids_col)

                for ctrl_type in [self.neg_ctrltype, 'pos', 'exp']:
                    df_2_labeled_ctrltype = df_2_labeled[df_2_labeled.ctrl_type==ctrl_type]
                    print(f"[{ctrl_type}] {col_suffix} hit rate : {round(df_2_labeled_ctrltype[label_col].mean(),4)} (n={len(df_2_labeled_ctrltype)})")
            else:
                df_2_labeled = None
                
        return df_1_labeled, df_2_labeled
    
    
    def get_binary_labels_from_lcms_data(self, df, split_by_list=[''], binary_label_col='measured_binary_score', nonbinary_label_col='prod_conc_lcms_actual'):
        
        df_labeled = df.copy()

        for col_suffix in split_by_list:
            nonbinary_label_col_wsuffix = f'{nonbinary_label_col}_{col_suffix}' if col_suffix!='' else nonbinary_label_col
            output_col_name = f'{binary_label_col}{self.col_suffix_db[col_suffix]}_lcms' if col_suffix!='' else f'{binary_label_col}_lcms'

            # get LCMS mean and stdev of negctrl samples -> get neg ctrl threshold
            df_neg = df.loc[df.ctrl_type=='EV', nonbinary_label_col_wsuffix]
            neg_mean = df_neg.median()
            neg_stdev = df_neg.quantile(0.75) - df_neg.quantile(0.5)
            thres_neg = (neg_mean + neg_stdev)*1.5

            # get LCMS mean and stdev of posctrl samples
            df_pos = df.loc[df.ctrl_type=='pos', nonbinary_label_col_wsuffix]
            pos_mean = df_pos.median()
            pos_stdev = df_pos.quantile(0.5) - df_pos.quantile(0.25)
            thres_pos = pos_mean - pos_stdev

            # get combined threshold and binarize nonbinary labels
            nonbinary_thres = max(thres_neg, thres_pos)
            print(f'[LCMS] neg ctrl thres: {round(thres_neg,4)}; pos ctrl thres: {round(thres_pos,4)}; FINAL THRES: {round(nonbinary_thres,4)}')
            binary_labels = (df[nonbinary_label_col_wsuffix] > nonbinary_thres)*1
            df_labeled.loc[:, output_col_name] = binary_labels

            print('EV control mean:', round(df_labeled.loc[df_labeled.ctrl_type=='EV', output_col_name].mean(),3))
            print('pos control mean:', round(df_labeled.loc[df_labeled.ctrl_type=='pos', output_col_name].mean(),3))
            print('exp sample mean:', round(df_labeled.loc[df_labeled.ctrl_type=='exp', output_col_name].mean(),3))
        return df_labeled
    

    def update_maldi_binary_labels_w_lcms_labels(self, df_maldi_lcms, split_by_list=[''], neg_ctrltype='EV',
                                          binary_label_col_base='true_score_binary', nonbinary_label_col_base='prod_conc_lcms_actual'):
        
        # get lcms-based binary labels
        df_maldi_lcms_labeled = self.get_binary_labels_from_lcms_data(df_maldi_lcms, split_by_list, binary_label_col=binary_label_col_base, nonbinary_label_col=nonbinary_label_col_base)
            
        for col_suffix in split_by_list:
            # get col names
            binary_label_col_lcms = f'{binary_label_col_base}{self.col_suffix_db[col_suffix]}' if col_suffix!='' else binary_label_col_base
            binary_label_col = f'{binary_label_col_base}_{col_suffix}' if col_suffix!='' else binary_label_col_base
            nonbinary_label_col = f'{nonbinary_label_col_base}{self.col_suffix_db[col_suffix]}' if col_suffix!='' else nonbinary_label_col_base

            # combine MALDI and LCMS binary labels (if LCMS or MALDI are 0, final label should be 0) and update maldi_lcms dataframe with new column
            df_maldi_lcms_binary_labels = df_maldi_lcms_labeled[f'{binary_label_col_lcms}_lcms'].to_numpy()
            df_maldi_lcms_binary_labels_old = df_maldi_lcms_labeled[binary_label_col].to_numpy()
            df_maldi_lcms_binary_labels_new = df_maldi_lcms_binary_labels_old * df_maldi_lcms_binary_labels
            df_maldi_lcms_labeled[f'{binary_label_col}_lcms'] = list(df_maldi_lcms_binary_labels_new)
        
        return df_maldi_lcms_labeled
    
    
    def propagate_labels_to_maldi_detections_2(self, df_maldi_lcms_labeled, df_2, split_by_list=[''], binary_label_col_base='true_score_binary'):
        
        for col_suffix in split_by_list:
            
            # get col names
            binary_label_col_lcms = f'{binary_label_col_base}{self.col_suffix_db[col_suffix]}' if col_suffix!='' else binary_label_col_base
            binary_label_col = f'{binary_label_col_base}_{col_suffix}' if col_suffix!='' else binary_label_col_base
            plate_col = f'plate_{col_suffix}' if col_suffix!='' else 'plate'
            address_col = f'address_{col_suffix}' if col_suffix!='' else 'address'

            # propagate these labels back to df_2 by merging by maldi plate and address
            df_maldi_lcms_labeled_tomerge = df_maldi_lcms_labeled[[plate_col, address_col, f'{binary_label_col_lcms}_lcms']]
            df_maldi_lcms_labeled_tomerge = df_maldi_lcms_labeled_tomerge[~df_maldi_lcms_labeled_tomerge[plate_col].isnull()]

            # replace <binary_label_col> with <binary_label_col>_lcms vals; store old values in <binary_label_col>_maldi
            df_2_labeled = df_2.copy()
            df_2_labeled = df_2.merge(df_maldi_lcms_labeled_tomerge, how='inner', on=[plate_col, address_col])
            df_maldi_lcms_labeled[f'{binary_label_col}_maldi'] = df_maldi_lcms_labeled[binary_label_col]
            df_maldi_lcms_labeled[binary_label_col] = df_maldi_lcms_labeled[f'{binary_label_col_lcms}_lcms']
            df_2_labeled[f'{binary_label_col}_maldi'] = df_2_labeled[binary_label_col]
            df_2_labeled[binary_label_col] = df_2_labeled[f'{binary_label_col_lcms}_lcms']
            mean_binary_label_old = round(df_2_labeled[f'{binary_label_col}_maldi'].mean(), 5)
            mean_binary_label_new = round(df_2_labeled[f'{binary_label_col_lcms}_lcms'].mean(), 5)
            print(f'[{col_suffix}] Average MALDI 2-sample dataset label: {mean_binary_label_old} (MALDI thresholds); {mean_binary_label_new} (MALDI + LCMS thresholds)')
        
        return df_2_labeled
    
    
    def propagate_labels_to_maldi_detections_1(self, df_2_labeled, df_1, split_by_list=[''], binary_label_col_base='true_score_binary'):
        
        for col_suffix in split_by_list:
            
            # get col names
            binary_label_col_lcms = f'{binary_label_col_base}{self.col_suffix_db[col_suffix]}' if col_suffix!='' else binary_label_col_base
            binary_label_col = f'{binary_label_col_base}_{col_suffix}' if col_suffix!='' else binary_label_col_base
            spectra_ids_col = f'spectra_ids_{col_suffix}' if col_suffix!='' else 'spectra_ids'
            # propagate labels from df_2 back to df_1
            df_1_labeled = df_1.copy()
            df_1_labeled[f'{binary_label_col}_lcms'] = np.nan
            
            for i in list(df_2_labeled.index):
                binary_label_lcms = df_2_labeled.loc[i, f'{binary_label_col_lcms}_lcms']
                spectra_ids = str(df_2_labeled.loc[i, spectra_ids_col]).split(', ')
                binary_labels_old = df_1_labeled.loc[df_1_labeled[spectra_ids_col].isin(spectra_ids), binary_label_col].to_numpy()
                binary_labels_new = binary_labels_old * binary_label_lcms
                df_1_labeled.loc[df_1_labeled[spectra_ids_col].isin(spectra_ids), f'{binary_label_col_lcms}_lcms'] = list(binary_labels_new)
            df_1_labeled[f'{binary_label_col}_maldi'] = df_1_labeled[binary_label_col]
            df_1_labeled[binary_label_col] = df_1_labeled[f'{binary_label_col_lcms}_lcms']
            mean_binary_label_old = round(df_1_labeled[f'{binary_label_col}_maldi'].mean(), 5)
            mean_binary_label_new = round(df_1_labeled[f'{binary_label_col_lcms}_lcms'].mean(), 5)
            print(f'[{col_suffix}] Average MALDI 1-sample dataset label: {mean_binary_label_old} (MALDI thresholds); {mean_binary_label_new} (MALDI + LCMS thresholds)')

        return df_1_labeled

    
    def dataframe_to_XYdataset(self, df, X_colnames, y_colname=None, col_suffix='', dropna=True):
        if col_suffix != '':
            col_suffix = f'_{col_suffix}'
        X_colnames = [f'{x_colname}{col_suffix}' for x_colname in X_colnames]
        
        if y_colname is not None:
            y_colname = f'{y_colname}{col_suffix}'
            df = df[X_colnames+[y_colname]]
            if dropna:
                df = df.dropna()
            X = np.array(df[X_colnames])
            y = np.array(list(df[y_colname]))
            if y[0] == None: y = None
        else:
            df = df[X_colnames]
            if dropna:
                df = df.dropna()
            X = np.array(df[X_colnames])
            y = None
        return X, y
    
        
    def train_maldi2maldi_binary(self, df, X_colnames, y_colname, col_suffix='', class_weight={0:1,1:1},
                                models={
                                    'LogisticRegression(L2)': LogisticRegression(),
                                    'LogisticRegression(L1)': LogisticRegression(penalty='l1',solver='liblinear'),
                                    'RandomForestClassifier':RandomForestClassifier(n_estimators=10),
                                    'AdaBoostClassifier': AdaBoostClassifier(n_estimators=200),
#                                     'MLPClassifier': MLPClassifier(hidden_layer_sizes=(50,25,12,4), activation='relu', learning_rate_init=1e-3, max_iter=2000, shuffle=True)
                                        }):
        """
        Train classification models
        """
        
        # get dataset
        X, y = self.dataframe_to_XYdataset(df, X_colnames, y_colname, col_suffix)
        unique_labels, counts = np.unique(y[:], return_counts=True)
        print(unique_labels, counts)
        print('Dataset shape:', X.shape, '; Average label:', round(np.mean(y), 4))
        
        if len(unique_labels)>1:
            if models is None:
                models = {
                            'LogisticRegression(L2)': LogisticRegression(class_weight=class_weight),
                            'LogisticRegression(L1)': LogisticRegression(penalty='l1',solver='liblinear', class_weight=class_weight),
                            'RandomForestClassifier':RandomForestClassifier(n_estimators=10),
                            'AdaBoostClassifier': AdaBoostClassifier(n_estimators=200),
#                             'MLPClassifier': MLPClassifier(hidden_layer_sizes=(50,25,12,4), activation='relu', learning_rate_init=1e-3, max_iter=2000, shuffle=True)
                }

            try:
                models, model_metrics_cv = get_ensemble_cv_results(X, y, models, n_splits=5, verbose=True, classifier_or_regressor=0)
            except:
                models, model_metrics_cv = get_ensemble_cv_results(X, y, models, n_splits=None, verbose=True, classifier_or_regressor=0)
        else:
            models, model_metrics_cv = None, None

        return models, model_metrics_cv
    
    
    def train_maldi2lcms_nonbinary(self, df, X_colnames, y_colname, col_suffix='',
                                   models={
                                    #'LinearRegression': LinearRegression(),
                                    # 'RidgeRegression': Ridge(),
                                      'LassoRegression': Lasso(),
                                      'RandomForestRegressor':RandomForestRegressor(n_estimators=10),
                                      'AdaBoostRegressor': AdaBoostRegressor(n_estimators=10),
                                      'MLPRegressor': MLPRegressor(hidden_layer_sizes=(50,25,12,4), activation='relu', learning_rate_init=1e-3, max_iter=10000, shuffle=True)
                                   }):
        """
        Train regression models
        """
        # get dataset
        X, y = self.dataframe_to_XYdataset(df, X_colnames, y_colname, col_suffix)
        print('Dataset shape:', X.shape, '; Average label:', round(np.mean(y), 4))
        
        # train models
        models, model_metrics_cv = get_ensemble_cv_results(X, y, models, n_splits=5, verbose=True, classifier_or_regressor=1)
        
        return models, model_metrics_cv
    
    
    def predict_maldi2maldi_binary(self, models, df, X_colnames, y_colname, col_suffix=''):
        """
        Predict using classification models
        """
        # get dataset
        X, y = self.dataframe_to_XYdataset(df, X_colnames, y_colname, col_suffix, dropna=False)
        print('Dataset shape:', X.shape)
        if y is not None: print('Average label:', round(np.nanmean(y), 4))
            
        if models is not None:
            # perform predictions
            nonan_idx = np.where(~np.isnan(X[:,0]))[0]
            X_nonan = X[nonan_idx,:]
            y_pred_nonan, y_pred_variance_nonan = get_ensemble_predictions(models, X_nonan,
                                                                            classifier_or_regressor=0, verbose=True)
            y_pred, y_pred_variance = np.empty((len(X),)), np.empty((len(X),))
            y_pred[:] = np.nan
            y_pred_variance[:] = np.nan
            y_pred[nonan_idx] = y_pred_nonan
            y_pred_variance[nonan_idx] = y_pred_variance_nonan
        else:
            y_pred, y_pred_variance = np.empty((len(X),)), np.empty((len(X),))
            y_pred[:] = np.nan
            y_pred_variance[:] = np.nan
            
        return y_pred, y_pred_variance
    
    
    def predict_maldi2lcms_nonbinary(self, models, df, X_colnames, y_colname, col_suffix=''):
        """
        Predict using classification models
        """
        # get dataset
        X, y = self.dataframe_to_XYdataset(df, X_colnames, y_colname, col_suffix, dropna=False)
        print('Dataset shape:', X.shape)
        if y is not None: print('Average label:', round(np.nanmean(y), 4))
        
        # perform predictions
        nonan_idx = np.where(~np.isnan(X[:,0]))[0]
        X_nonan = X[nonan_idx,:]
        y_pred_nonan, y_pred_variance_nonan = get_ensemble_predictions(models, X_nonan,
                                                                             classifier_or_regressor=1, verbose=True)
        y_pred, y_pred_variance = np.empty((len(X),)), np.empty((len(X),))
        y_pred[:] = np.nan
        y_pred_variance[:] = np.nan
        y_pred[nonan_idx] = y_pred_nonan
        y_pred_variance[nonan_idx] = y_pred_variance_nonan
        
        return y_pred, y_pred_variance

    
    def select_variants(self, selected_variants,
                        select_by=['predicted_binary_score', 'predicted_nonbinary_score', 'predicted_nonbinary_score_(r)'],
                        max_num_hits=100, select_up_to_max=False, col_suffix_list=['']):
        """
        Select hits by ranking predicted LCMS scores and MALDI hits
        """
        
        # perform hits selection based on metric
        selected_variants_all = selected_variants.copy()
        for metric in select_by:

            # binary score selection
            if metric in ['predicted_binary_score', 'predicted_nonbinary_score', 'measured_binary_score']:
                hit_selection_conditions = {f"selected_variants_all['{metric}{self.col_suffix_db[col_suffix]}']": (0, '>') for col_suffix in col_suffix_list}
                selected_variants_all = selected_variants_all[eval(get_metric_filter_str(hit_selection_conditions, joiner=' | '))]
                print(f'Selected {len(selected_variants_all)} variants based on {metric}.')

            # pick top {max_num_hits} hits from nonbinary measure
            elif metric in ['enantiomeric_excess_(plus-minus)', 'predicted_nonbinary_score_(r)', 'predicted_nonbinary_score_(+)', 'predicted_nonbinary_score_(-)', 'measured_nonbinary_score_(r)', 'measured_nonbinary_score_(+)', 'measured_nonbinary_score_(-)']:
                selected_variants_all = selected_variants_all.sort_values(by=[metric], ascending=False)
                if max_num_hits is not None and len(selected_variants_all)>max_num_hits:
                    selected_variants_all = selected_variants_all.iloc[:max_num_hits]
                    print(f'Selected {len(selected_variants_all)} variants based on {metric}.')
            
        print(f'# of variants selected based on metrics filtering: {len(selected_variants_all)}')
        
        # determine whether to select top <n> variants (regardless of hit status), or just the hits
        if select_up_to_max:
            max_num_hits = min(max_num_hits, len(selected_variants))
            selected_variants_all = selected_variants.copy().sort_values(by=[f'{metric}{self.col_suffix_db[col_suffix]}' for col_suffix in col_suffix_list], ascending=False)
            selected_variants_all = selected_variants_all[:max_num_hits]
            
        print(f'FINAL # of variants selected: {len(selected_variants_all)}')

        return selected_variants_all
    
    
    def generate_worklist(self, selected_variants, track_dilution_plates=False, col_suffix='', plate_type_to_select_from='maldi_plate'):
        
        from utils.additional_utils import create_alphabet_integer_mapping, address_to_xy
        from utils.table_properties import lcms_worklist_columns

        col_suffix_db = self.col_suffix_db[col_suffix]
        
        # get MALDI plate barcode
        plate = list(selected_variants[f'{plate_type_to_select_from}{col_suffix_db}'])
        plate = [plate_str.split(',')[0] for plate_str in plate][0]
        
        # get unique enzyme barcodes
        enzyme_barcodes = list(selected_variants.enzyme_barcode)
        
        # get unique source plate barcodes
        source_wells = list(selected_variants[f'source_address{col_suffix_db}'])
        source_plates = list(selected_variants[f'source_plate{col_suffix_db}'])
        unique_source_plates = list(set(source_plates))
        lcms_final_plates = [self.map_maldi_to_lcms_plate(source_plate)['lcms'] for source_plate in unique_source_plates]
        source_plate_to_final_plate_mapping = {sp:fp for sp, fp in zip(unique_source_plates, lcms_final_plates)}
        print('LCMS intermediate to final plate mapping: ', source_plate_to_final_plate_mapping)
        
        # initialize worklist
        worklist = pd.DataFrame(index=list(range(len(enzyme_barcodes))), columns=lcms_worklist_columns)
        worklist.index.name = 'Sample ID'

        # Standard entries
        worklist.loc[:, 'Sample Name'] = list(selected_variants[f'substrate_barcode{col_suffix_db}'])
        worklist['Plate Code'] = 'PlateOrVial'
        worklist['Sample Type'] = 'Sample'
        worklist['Level Name'] = ''
        worklist['Inj Vol'] = 'As Method'
        worklist['Dilution'] = 1
        worklist['Comment'] = ''
        worklist['Info.'] = ''
        worklist['Sample Group'] = ''

        integer_to_alphabet_dict = dict((v,k) for k,v in create_alphabet_integer_mapping().items())
        for i in list(worklist.index):
            
            source_plate, source_well = source_plates[i], source_wells[i]
            
            if source_well is not None:
                source_well = source_well.split(',')[0] # take only the first well if there are multiple
                x, y = address_to_xy(source_well, zero_index=False)
                well_position = f"{integer_to_alphabet_dict[y]}{x}"
            else:
                well_position = enzyme_barcodes[i]
                # add code to get wel position from enzyme barcode

            # Barcode (get final dilution plate barcode from rxn source plate using transfer tracking)
            if source_plate is not None and plate_type_to_select_from=='maldi_plate':
                lcms_plate = source_plate_to_final_plate_mapping[source_plate]
            else:
                lcms_plate = plate

            # Sample Position
            sample_position = well_position

            # Method
            method = ''

            # Data File
            data_file = f"{str(datetime.now())[:11].replace('-','')}_{lcms_plate}_{well_position}.d"

            # update worklist
            worklist.loc[i, ['Barcode', 'Sample Position', 'Method', 'Data File']] = [lcms_plate, sample_position, method, data_file]

        # post worklist
        r = s3_df2csv(worklist, 'lcms-file-store', f'worklists/{plate}{col_suffix_db}_LCMSworklist.csv')
        print(f"'{plate}{col_suffix_db}_LCMSworklist.csv' posted to S3 lcms-file-store bucket.")
        
        return worklist
    
    
    def map_maldi_to_lcms_plate(self, source_plate):
        from utils.lims_utils import track_transfers_in_out, LIMS_BASE_URL_PROD, LIMS_BASE_URL_DEV
        plate_sequence = {'source':source_plate}
        
        try:
            # get LCMS final plate
            _, transfers_out = track_transfers_in_out(source_plate, well_details=None, in_strings=None, out_strings=['lcms'], lims_base_url=LIMS_BASE_URL_PROD)
            if 'X12Y8' in transfers_out:
                lcms_plate_barcode = transfers_out['X12Y8']['plate']
                plate_sequence.update({'lcms':lcms_plate_barcode[0]})
            else:
                print(f'No LCMS final plates found for source plate {source_plate}.')
                self.message_suffix += f'No LCMS final plates found for source plate {source_plate}. '
                plate_sequence.update({'lcms':source_plate})
        except:
            print(f'No LCMS final plates found for source plate {source_plate}.')
            self.message_suffix += f'No LCMS final plates found for source plate {source_plate}.'
            plate_sequence.update({'lcms':source_plate})
        return plate_sequence
    
    
    def consolidate_table_columns(self, df, cols_to_consolidate, suffix=''):
        cols_to_consolidate_wsuffix = [f'{c}{suffix}' for c in cols_to_consolidate if f'{c}{suffix}' in df.columns.tolist()]
        # combine content in columns
        fill_cols = df.loc[:, cols_to_consolidate_wsuffix].rename(columns={col_wsuffix: col for col,col_wsuffix in zip(cols_to_consolidate, cols_to_consolidate_wsuffix)})
        cols_tofill = df.loc[:, cols_to_consolidate]
        cols_tofill = cols_tofill.fillna(fill_cols)
        df.loc[:, cols_to_consolidate] = cols_tofill
        # remove duplicate columns
        df = df.drop(columns=cols_to_consolidate_wsuffix)
        return df

    def GET_MALDI_DATASETS(self,
                           plate=None, run=None, exp_workflow=None, maldi_1_csv=None, maldi_2_csv=None,
                           combine_by=['address', 'source_plate_grp'], split_by_list=None, get_split_by_list=False,
                           source_plate_groupings={}, get_source_plate_groupings=False, plate_to_grp_idx_mapping={},
                           data_types_to_get=['maldi_detections_1','maldi_detections_2'], maldi_cols_to_retain=maldi_basic_metadata, model_type_list=[],
                           get_binary_labels_maldi=False):
        
        self.source_plate_groupings = source_plate_groupings
        self.plate_to_grp_idx_mapping = plate_to_grp_idx_mapping
        
        ##########################################
        # get 1-sample + 2-sample MALDI datasets #
        ##########################################
        
        data = {}
        if ('maldi1sample_classification' in model_type_list) or ('maldi_detections_1' in data_types_to_get):
            data[self.maldi_detections_1] = self.get_maldi_1sample_data(plate, run, exp_workflow, maldi_1_csv,
                                                                        get_source_plate_groupings=get_source_plate_groupings,
                                                                        source_plate_groupings=self.source_plate_groupings,
                                                                        plate_to_grp_idx_mapping=self.plate_to_grp_idx_mapping)
            self.plate_list = list(set(data[self.maldi_detections_1].plate))
        if ('maldi2sample_classification' in model_type_list) or ('maldi2lcms_regression' in model_type_list) or ('maldi_detections_2' in data_types_to_get):
            data[self.maldi_detections_2] = self.get_maldi_2sample_data(plate, run, exp_workflow, maldi_2_csv, data[self.maldi_detections_1],
                                                                        get_source_plate_groupings=get_source_plate_groupings,
                                                                        source_plate_groupings=self.source_plate_groupings,
                                                                        plate_to_grp_idx_mapping=self.plate_to_grp_idx_mapping)
            self.plate_list = list(set(data[self.maldi_detections_2].plate))
            
        ############################################################
        # get substrates in dataset & filter dataset by substrates #
        ############################################################
        
        if get_split_by_list:
            if split_by_list is None:
                self.split_by_list = ['']
            else:
                if split_by_list==[]:
                    substrate_list = []
                    for k, v in self.substrate_list.items():
                        substrate_list += v
                    substrate_list = sorted(list(set(substrate_list)))
                    self.split_by_list = substrate_list.copy()
                else:
                    self.split_by_list = split_by_list
                    
            split_by_list = self.split_by_list
        
        ###############################################################
        # merge MALDI data based on common address & source plate grp #
        ###############################################################
        
        if split_by_list != ['']:
            for data_type in [self.maldi_detections_1, self.maldi_detections_2]:
                # filter dataset by substrate
                data[data_type] = data[data_type][data[data_type].substrate_barcode.isin(split_by_list)]
                df = data[data_type]
                substrate_list_datatype = sorted(list(set(df.substrate_barcode)))
                print(data_type, f'{len(data[data_type])} rows. {len(substrate_list_datatype)} substrates: {substrate_list_datatype}')
                data[data_type] = combine_rows_by_unique_enzyme(df, cols_to_retain=maldi_cols_to_retain,
                                                                colnames_to_modify=maldi_additional_metadata + maldi_detections_derived,
                                                                combine_by=combine_by,
                                                                col_to_split_by='substrate_barcode',
                                                                split_by_list=split_by_list,
                                                                test_colname='ptr')
                if get_source_plate_groupings:
                    data[data_type] = data[data_type].astype({"source_plate_grp": int})
                print(data_type, f'{len(data[data_type])} rows.')
            self.substrate_list_maldi = substrate_list_datatype

        ####################
        # label maldi data #
        ####################
        if get_binary_labels_maldi:
            data[self.maldi_detections_1], data[self.maldi_detections_2] = self.get_binary_labels_from_maldi_data(
                data[self.maldi_detections_1], data[self.maldi_detections_2],
                split_by_list=split_by_list, neg_ctrltype=self.neg_ctrltype, maldi_label_condition=self.maldi_label_condition
                )
            print(f'Updated MALDI dataset labels.')
            
        ###############################
        ## save MALDI datasets in S3 ##
        ###############################
        for data_type, df in data.items():
            s3_df2csv(df, self.s3_bucket, f'{self.s3_subfolder}{data_type}{self.fname_suffix}.csv')
        
        return data
    
    
    def GET_LCMS_DATASETS(self, lcms_run=None, exp_workflow=None, lcms_csv=None, lcms_chiral_csv=None, combine_by=['address', 'source_plate_grp'],
                              source_plate_groupings={}, get_source_plate_groupings=True, plate_to_grp_idx_mapping={}, split_by_list=None,
                              data_types_to_get=['lcms_detections','lcms_detections_chiral'],
                              C18_cols_to_retain=lcms_basic_metadata + library_metadata + plate_time_metadata + pelletqc_metadata,
                              chiral_merge_method='outer', calculate_conversion_enantioselectivity=True,
                              perform_qc=False, plate_col='plate', subconc_init=100,
                              get_binary_labels_lcms=False, nonbinary_label_col='prod_conc_lcms_actual'):
        
        self.source_plate_groupings = source_plate_groupings
        self.plate_to_grp_idx_mapping = plate_to_grp_idx_mapping
        data = {}
        
        ####################################
        # get LCMS C18 and chiral datasets #
        ####################################
        
        if exp_workflow is None:
            if isinstance(lcms_run, list):
                lcms_run_C18 = [f'{r}_C18' for r in lcms_run]
                lcms_run_chiral = [f'{r}_chiral' for r in lcms_run]
            else:
                lcms_run_C18 = f'{lcms_run}_C18'
                lcms_run_chiral = f'{lcms_run}_chiral'
        else:
            lcms_run_C18 = None
            lcms_run_chiral = None

        if 'lcms_detections' in data_types_to_get:
            data[self.lcms_detections] = self.get_lcms_C18_data(lcms_run_C18, exp_workflow, lcms_csv,
                                                                get_source_plate_groupings=get_source_plate_groupings,
                                                                source_plate_groupings=self.source_plate_groupings,
                                                                plate_to_grp_idx_mapping=self.plate_to_grp_idx_mapping)
        if 'lcms_detections_chiral' in data_types_to_get:
            data[self.lcms_detections_chiral] = self.get_lcms_chiral_data(lcms_run_chiral, exp_workflow, lcms_chiral_csv,
                                                                get_source_plate_groupings=get_source_plate_groupings,
                                                                source_plate_groupings=self.source_plate_groupings,
                                                                plate_to_grp_idx_mapping=self.plate_to_grp_idx_mapping)
        
        ##################################################################
        # merge LCMS C18 data based on common address & source plate grp #
        ##################################################################
        
        if split_by_list != ['']:
            # filter dataset by substrate
            data_type = self.lcms_detections
            data[data_type] = data[data_type][data[data_type].substrate_barcode.isin(split_by_list)]
            df = data[data_type]
            self.substrate_list_lcms = sorted(list(set(df.substrate_barcode)))
            print(data_type, f'{len(data[data_type])} rows. {len(self.substrate_list_lcms)} substrates: {self.substrate_list_lcms}')
            data[data_type] = combine_rows_by_unique_enzyme(df, cols_to_retain=C18_cols_to_retain,
                                                    colnames_to_modify=lcms_detections_area + lcms_additional_metadata + lcms_detections_derived,
                                                            combine_by=combine_by,
                                                            col_to_split_by='substrate_barcode',
                                                            split_by_list=self.substrate_list_lcms,
                                                            test_colname='prod_conc_lcms_actual')
            if get_source_plate_groupings:
                print(data[data_type].columns.tolist())
                data[data_type] = data[data_type].astype({"source_plate_grp": int})
            print(data_type, f'{len(data[data_type])} rows.')
        
        ########################################
        ## combine lcms C18 & chiral datasets ##
        ########################################
        
        self.data = data

        if self.lcms_detections in data and self.lcms_detections_chiral in data:
            print(self.lcms_detections_chiral, f'{len(data[self.lcms_detections_chiral])} rows.')
            
            # get columns to combine by -- modify lcms C18 table to contain these columns
            for c in self.combine_maldi_lcms_by:
                if c not in data[self.lcms_detections]:
                    cmod = f'{c}_{self.split_by_list[0]}' if self.split_by_list[0] != '' else c
                    data[self.lcms_detections].loc[:,c] = data[self.lcms_detections][cmod].tolist()
                    
            # get lcms chiral columns to include
            lcms_cols_chiral = [] 
            for c in self.combine_maldi_lcms_by + common_metadata + lcms_detections_chiral_area + lcms_detections_chiral_derived + library_metadata + plate_time_metadata + pelletqc_metadata:
                if c not in lcms_cols_chiral:
                    lcms_cols_chiral.append(c)
            
            # perform merge
            data['lcms_detections_all'] = data[self.lcms_detections].merge(
                data[self.lcms_detections_chiral][lcms_cols_chiral],
                how=chiral_merge_method,
                on=self.combine_maldi_lcms_by,
                suffixes=['','_chiral']
            )
            print(f'Added LCMS chiral data to LCMS C18 dataset using columns {self.combine_maldi_lcms_by}. {len(data["lcms_detections_all"])} rows.')
            
            # consolidate columns
            data['lcms_detections_all'] = self.consolidate_table_columns(data['lcms_detections_all'].copy(), self.cols_to_consolidate['C18_chiral'], suffix='_chiral')
            data['lcms_detections_all'] = data['lcms_detections_all'].drop(columns=['id', 'address_chiral', 'plate'])
            
            # remove rows where 'source_plate_group' or 'source_address' are null
            for c in self.combine_maldi_lcms_by:
                data['lcms_detections_all'] = data['lcms_detections_all'][~data['lcms_detections_all'][c].isnull()]
                
            # cast hamming column as integer type
            data['lcms_detections_all']['hamming'] = pd.to_numeric(data['lcms_detections_all']['hamming'], downcast='integer')
                        
        elif self.lcms_detections in data and self.lcms_detections_chiral not in data:
            data['lcms_detections_all'] = data[self.lcms_detections].copy()
            
        elif self.lcms_detections not in data and self.lcms_detections_chiral in data:
            data['lcms_detections_all'] = data[self.lcms_detections_chiral].copy()
            
           
        #################################################
        ## get derived LCMS values for ee & conversion ##
        #################################################
        
        ## C18 calculations ##
        if self.lcms_detections in data and calculate_conversion_enantioselectivity:
            conversion, ee_plusoverminus, eeratio_plusoverminus = calculate_conversion_enantioselectivity_C18(
                data['lcms_detections_all'], split_by_list, enantiomer_barcodes=self.plot_enantiomer_lcms_metrics
            )
            for s in conversion:
                data['lcms_detections_all'].loc[:,f'measured_conversion{self.col_suffix_db[s]}_C18'] = conversion[s]
            data['lcms_detections_all'].loc[:,'measured_enantiomeric_excess_(+over-)_C18'] = ee_plusoverminus
            data['lcms_detections_all'].loc[:,'measured_enantiomeric_ratio_(+over-)_C18'] = eeratio_plusoverminus
            print(f'Calculated conversion and enantioselectivity metrics for C18 LCMS data. ')
        
        ## Chiral calculations ##
        if self.lcms_detections_chiral in data and calculate_conversion_enantioselectivity:
            conversion, ee_plusoverminus, eeratio_plusoverminus = calculate_conversion_enantioselectivity_chiral(data['lcms_detections_all'], ee_calc_method=3)
            for s in conversion:
                data['lcms_detections_all'].loc[:,f'measured_conversion{s}_chiral'] = conversion[s]
            data['lcms_detections_all'].loc[:,'measured_enantiomeric_excess_(+over-)_chiral'] = ee_plusoverminus
            data['lcms_detections_all'].loc[:,'measured_enantiomeric_ratio_(+over-)_chiral'] = eeratio_plusoverminus
            print(f'Calculated conversion and enantioselectivity metrics for chiral LCMS data. ')
        
        ################################################
        ## update derived LCMS derived metric columns ##
        ################################################
        if self.lcms_detections in data or self.lcms_detections_chiral in data_types_to_get:

            # MEASURED_NONBINARY_SCORE & CONVERSION_(R)
            for col_suffix in self.substrate_list_lcms:
                data['lcms_detections_all'][f'measured_nonbinary_score_{col_suffix}'] = data['lcms_detections_all'][f'{self.lcms_output_feature}_{col_suffix}'].tolist()
            data['lcms_detections_all'].loc[:, 'measured_conversion_(r)'] = data['lcms_detections_all'][f'measured_conversion_(r)_C18'].tolist()
                
            # MEASURED CONVERSION & ENANTIOSELECTIVITY METRICS
            if self.lcms_detections_chiral in data:
                method_suffix = '_chiral'
            elif self.lcms_detections in data:
                method_suffix = '_C18'
            if calculate_conversion_enantioselectivity:
                if f'measured_conversion_(+){method_suffix}' in data['lcms_detections_all']:
                    data['lcms_detections_all'].loc[:, 'measured_enantiomeric_ratio_(+over-)'] = data['lcms_detections_all'][f'measured_enantiomeric_ratio_(+over-){method_suffix}'].tolist()
                    data['lcms_detections_all'].loc[:, 'measured_enantiomeric_excess_(+over-)'] = data['lcms_detections_all'][f'measured_enantiomeric_excess_(+over-){method_suffix}'].tolist()
                    data['lcms_detections_all'].loc[:, 'measured_conversion_(+)'] = data['lcms_detections_all'][f'measured_conversion_(+){method_suffix}'].tolist()
                    data['lcms_detections_all'].loc[:, 'measured_conversion_(-)'] = data['lcms_detections_all'][f'measured_conversion_(-){method_suffix}'].tolist()
                else:
                    data['lcms_detections_all'].loc[:, 'measured_enantiomeric_ratio_(+over-)'] = np.nan
                    data['lcms_detections_all'].loc[:, 'measured_enantiomeric_excess_(+over-)'] = np.nan
                    data['lcms_detections_all'].loc[:, 'measured_conversion_(+)'] = np.nan
                    data['lcms_detections_all'].loc[:, 'measured_conversion_(-)'] = np.nan
    
        ######################################
        ## get binary labels from lcms data ##
        ######################################
        if get_binary_labels_lcms:
            data['lcms_detections_all'] = self.get_binary_labels_from_lcms_data(data['lcms_detections_all'].copy(), nonbinary_label_col=nonbinary_label_col, split_by_list=split_by_list)
            print('Obtained binary labels from LCMS data.')
                
        if 'lcms_detections_chiral' in data_types_to_get:

            ######################
            # get lcms qc labels #
            ######################
            if perform_qc:
                data['lcms_detections_all'] = get_lcms_qc_derivedmetrics_checks(data['lcms_detections_all'], neg_ctrltype=self.neg_ctrltype, plate_col=plate_col, subconc_init=subconc_init)

            ##############################
            ## save LCMS datasets in S3 ##
            ##############################
            for data_type, df in data.items():
                s3_df2csv(df, self.s3_bucket, f'{self.s3_subfolder}{data_type}{self.fname_suffix}.csv')

        return data
    
    
    def GET_COMBINED_DATASETS(self,
                              plate=None, run=None, exp_workflow=None,
                              maldi_1_csv=None, maldi_2_csv=None, lcms_csv=None, lcms_chiral_csv=None,
                              source_plate_groupings=None, get_source_plate_groupings=False,
                              split_by_list=None, plate_to_grp_idx_mapping=None, get_split_by_list=False,
                              data_types_to_get=['maldi_detections_1','maldi_detections_2','lcms_detections','lcms_detections_chiral'], model_type_list=[],
                              chiral_merge_method='outer', get_binary_labels_maldi=False, get_binary_labels_lcms=False, calculate_conversion_enantioselectivity=False
                             ):

        self.source_plate_groupings = source_plate_groupings
        self.plate_to_grp_idx_mapping = plate_to_grp_idx_mapping
        
        data = {}
        
        #####################
        # get MALDI dataset #
        #####################
        data_maldi = self.GET_MALDI_DATASETS(plate=plate, run=run, exp_workflow=exp_workflow, maldi_1_csv=maldi_1_csv, maldi_2_csv=maldi_2_csv,
                                             split_by_list=split_by_list, get_split_by_list=get_split_by_list,
                                             get_source_plate_groupings=get_source_plate_groupings,
                                             source_plate_groupings=self.source_plate_groupings,
                                             plate_to_grp_idx_mapping=self.plate_to_grp_idx_mapping,
                                             data_types_to_get=data_types_to_get, model_type_list=model_type_list,
                                             get_binary_labels_maldi=get_binary_labels_maldi)
        data.update(data_maldi)
        
        ###################################
        # get LCMS C18 and chiral dataset #
        ###################################
        
        if 'lcms_detections' in data_types_to_get or 'lcms_detections_chiral' in data_types_to_get:
            data_lcms = self.GET_LCMS_DATASETS(lcms_run=run, exp_workflow=exp_workflow, lcms_csv=lcms_csv, lcms_chiral_csv=lcms_chiral_csv,
                                              get_source_plate_groupings=get_source_plate_groupings,
                                               source_plate_groupings=self.source_plate_groupings,
                                               plate_to_grp_idx_mapping=self.plate_to_grp_idx_mapping,
                                              split_by_list=split_by_list, data_types_to_get=data_types_to_get, chiral_merge_method=chiral_merge_method,
                                              calculate_conversion_enantioselectivity=calculate_conversion_enantioselectivity)

            data.update(data_lcms)
            
        #####################################
        ## combine maldi and lcms datasets ##
        #####################################
        
        if self.maldi_detections_2 in data and 'lcms_detections_all' in data:
            data[self.maldi_detections_2].astype({"source_plate_grp": int})
            data['lcms_detections_all'].astype({"source_plate_grp": int})
        
            for c in self.combine_maldi_lcms_by:
                if c not in data[self.maldi_detections_2]:
                    cmod = f'{c}_{self.split_by_list[0]}' if self.split_by_list[0] != '' else c
                    data[self.maldi_detections_2].loc[:,c] = data[self.maldi_detections_2][cmod].tolist()
            lcms_cols = data['lcms_detections_all'].columns.tolist()

            # perform merge
            data['maldi_lcms'] = data[self.maldi_detections_2].merge(
                data['lcms_detections_all'][lcms_cols],
                how='outer',
                on=self.combine_maldi_lcms_by,
                suffixes=['','_lcms']
            )
            # consolidate columns
            data['maldi_lcms'] = self.consolidate_table_columns(data['maldi_lcms'].copy(), self.cols_to_consolidate['maldi_lcms'], suffix='_lcms')
#             print(data['maldi_lcms'].columns.tolist())
            print(f'Combined 2-sample MALDI and LCMS datasets using columns {self.combine_maldi_lcms_by}.')

        else:
            data['maldi_lcms'] = data[self.maldi_detections_2].copy()

            ###############################################
            # label maldi_lcms data using lcms conditions #
            ###############################################

            if get_binary_labels_lcms:
                df_1_labeled = data[self.maldi_detections_1]
                df_2_labeled = data[self.maldi_detections_2] if self.maldi_detections_2 in data else None
                df_maldi_lcms_labeled = data['maldi_lcms']
                df_maldi_lcms_labeled = self.update_maldi_binary_labels_w_lcms_labels(df_maldi_lcms_labeled, split_by_list=split_by_list, binary_label_col_base='true_score_binary',
                                                                                      neg_ctrltype=self.neg_ctrltype, nonbinary_label_col_base=self.lcms_output_feature)
                if df_2_labeled is not None:
                    df_2_labeled = self.propagate_labels_to_maldi_detections_2(df_maldi_lcms_labeled, df_2_labeled, split_by_list=split_by_list,
                                                                               binary_label_col_base='true_score_binary')
                    df_1_labeled = self.propagate_labels_to_maldi_detections_1(df_2_labeled, df_1_labeled, split_by_list=split_by_list,
                                                                               binary_label_col_base='true_score_binary')
                else:
                    df_1_labeled = self.propagate_labels_to_maldi_detections_1(df_maldi_lcms_labeled, df_1_labeled, split_by_list=split_by_list,
                                                                               binary_label_col_base='true_score_binary')
                
                data[self.maldi_detections_1], data[self.maldi_detections_2], data['maldi_lcms'] = df_1_labeled, df_2_labeled, df_maldi_lcms_labeled
#             print(f'Updated MALDI & LCMS dataset labels for columns with suffix={col_suffix}.')
            
        ###################################
        ## save maldi_lcms dataset in S3 ##
        ###################################
        s3_df2csv(data['maldi_lcms'], self.s3_bucket, f'{self.s3_subfolder}maldi_lcms{self.fname_suffix}.csv')
            
        return data
    
    
    def GET_ENZYME_DATASET(self, data, split_by_list, predict_labels=False):
        
        # initialize
        enzyme_analytics_df = data['maldi_lcms'][['exp_workflow_barcode', 'exp_workflow_name', 'proj_barcode', 'proj_name', 'ctrl_type', 'enzyme_barcode', 'sequence', 'mutations', 'hamming','reference_enzyme']]

        # update metadata
        for col_suffix in split_by_list:
            col_suffix_db = self.col_suffix_db[col_suffix]
            col_suffix = col_suffix if col_suffix=='' else f'_{col_suffix}'
            enzyme_analytics_df.loc[:, f'maldi_run'] = list(data['maldi_lcms'][f'run{col_suffix}'])
            enzyme_analytics_df.loc[:, f'maldi_plate{col_suffix_db}'] = list(data['maldi_lcms'][f'plate{col_suffix}'])
            enzyme_analytics_df.loc[:, f'maldi_address{col_suffix_db}'] = list(data['maldi_lcms'][f'address{col_suffix}'])
            enzyme_analytics_df.loc[:, f'source_plate{col_suffix_db}'] = list(data['maldi_lcms'][f'source_plate{col_suffix}'])
            enzyme_analytics_df.loc[:, f'source_address{col_suffix_db}'] = list(data['maldi_lcms'][f'source_address{col_suffix}'])
            enzyme_analytics_df.loc[:, f'spectra_ids{col_suffix_db}'] = list(data['maldi_lcms'][f'spectra_ids{col_suffix}'])
            enzyme_analytics_df.loc[:, f'substrate_barcode{col_suffix_db}'] = list(data['maldi_lcms'][f'substrate_barcode{col_suffix}'])
            enzyme_analytics_df.loc[:, f'substrate_smiles{col_suffix_db}'] = list(data['maldi_lcms'][f'substrate_smiles{col_suffix}'])
            #### ADD LCMS PLATE METADATA UPDATE!!
        
        # update measured binary scores
        for col_suffix in self.substrate_list_maldi:
            col_suffix_db = self.col_suffix_db[col_suffix]
            enzyme_analytics_df.loc[:, f'measured_binary_score{col_suffix_db}'] = list(data['maldi_lcms'][f'true_score_binary_{col_suffix}'])
            
        # update measured nonbinary scores & sum
        if self.lcms_detections in data:
            for col_suffix in self.substrate_list_lcms:
                col_suffix_db = self.col_suffix_db[col_suffix]
                # update actual nonbinary scores
                enzyme_analytics_df.loc[:,f'measured_nonbinary_score{col_suffix_db}'] = list(data['maldi_lcms'][f'measured_nonbinary_score_{col_suffix}'])
                # update actual nonbinary sums
                if col_suffix_db=='_(r)':
                    enzyme_analytics_df.loc[:,f'measured_nonbinary_sum{col_suffix_db}'] = list(data['maldi_lcms'][f'sum_conc_lcms_actual_{col_suffix}'])
        
        # combine hits prediction from 1-sample and 2-sample classification to get final labels
        if predict_labels:
            for col_suffix in split_by_list:
                col_suffix_db = self.col_suffix_db[col_suffix]
                spectra_col = f'spectra_ids_{col_suffix}' if col_suffix!='' else 'spectra_ids'
                predicted_binary_score_col = f'predicted_binary_score_{col_suffix}' if col_suffix!='' else 'predicted_binary_score'
                predicted_binary_score_variance_col = f'predicted_binary_score_variance_{col_suffix}' if col_suffix!='' else 'predicted_binary_score_variance'
                predicted_nonbinary_score_col = f'predicted_nonbinary_score_{col_suffix}' if col_suffix!='' else 'predicted_nonbinary_score'
                predicted_nonbinary_score_variance_col = f'predicted_nonbinary_score_variance_{col_suffix}' if col_suffix!='' else 'predicted_nonbinary_score_variance'
                
                # obtain averaged binary labels from 1-sample and 2-sample predictions
                df_2_labeled = self.label_maldi2_with_maldi1(data['maldi_lcms'], data[self.maldi_detections_1],
                                                             predicted_binary_score_col, spectra_col)
                df_2_labeled = self.label_maldi2_with_maldi1(data['maldi_lcms'], data[self.maldi_detections_1],
                                                             predicted_binary_score_variance_col, spectra_col, binarize=False)
                labels_1sample = df_2_labeled[predicted_binary_score_col].to_numpy()
                labels_2sample = data['maldi_lcms'][predicted_binary_score_col].to_numpy()
                labels_1sample_variance = df_2_labeled[predicted_binary_score_variance_col].to_numpy()
                labels_2sample_variance = data['maldi_lcms'][predicted_binary_score_variance_col].to_numpy()

                # update predicted binary scores
                enzyme_analytics_df.loc[:, f'predicted_binary_score{col_suffix_db}'] = list(((labels_1sample*2 + labels_2sample)/3>=0.5)*1)
                enzyme_analytics_df.loc[:, f'predicted_binary_score_variance{col_suffix_db}'] = list((labels_1sample_variance*2 + labels_2sample_variance)/3)

                # update predicted nonbinary scores
                enzyme_analytics_df.loc[:, f'predicted_nonbinary_score{col_suffix_db}'] = list(data['maldi_lcms'][predicted_nonbinary_score_col])
                enzyme_analytics_df.loc[:, f'predicted_nonbinary_score_variance{col_suffix_db}'] = list(data['maldi_lcms'][predicted_nonbinary_score_variance_col])
                
        # calculate MEASURED enantiomeric excess -- applicable if non-racemic C18 LCMS or chiral LCMS
        if 'measured_nonbinary_score_(+)' in enzyme_analytics_df and self.plot_enantiomer_lcms_metrics is not None:
            enzyme_analytics_df.loc[:, f'measured_enantiomeric_excess_(+over-)'] = list(data['maldi_lcms'][f'measured_enantiomeric_excess_(+over-)'])
            enzyme_analytics_df.loc[:, f'measured_enantiomeric_ratio_(+over-)'] = list(data['maldi_lcms'][f'measured_enantiomeric_ratio_(+over-)'])
        
        # calculate PREDICTED enantiomeric excess -- applicable if non-racemic MALDI
        if len(split_by_list)>1 and self.plot_enantiomer_lcms_metrics is not None and predict_labels:
            nonbinary_score_plus = np.array(enzyme_analytics_df['predicted_nonbinary_score_(+)'])
            nonbinary_score_minus = np.array(enzyme_analytics_df['predicted_nonbinary_score_(-)'])
            enzyme_analytics_df.loc[:,f'predicted_enantiomeric_excess_(+over-)'] = list((nonbinary_score_plus - nonbinary_score_minus)/(nonbinary_score_plus + nonbinary_score_minus))
            enzyme_analytics_df.loc[:,f'predicted_enantiomeric_ratio_(+over-)'] = list(nonbinary_score_plus / nonbinary_score_minus)
        else:
            enzyme_analytics_df.loc[:,f'predicted_enantiomeric_excess_(+over-)'] = np.nan
            enzyme_analytics_df.loc[:,f'predicted_enantiomeric_ratio_(+over-)'] = np.nan
                
        return enzyme_analytics_df
    
    
    def TRAIN(self, plate=None, run=None, exp_workflow=None, lcms_method='',
              data_types_to_get=['maldi_detections_1', 'maldi_detections_2', 'lcms_detections', 'lcms_detections_chiral'],
              maldi_1_csv=None, maldi_2_csv=None, lcms_csv=None, lcms_chiral_csv=None, get_maldi_labels=False, get_lcms_labels=False,
              simulate_2sample_dataset=False, num_ref=16, num_rxn=4, num_unique_samples=48,
              get_source_plate_groupings=False, source_plate_groupings={}, plate_to_grp_idx_mapping={}, split_by_list=None,
              generate_heatmaps=True, skip_model_training=False):
        """
        Train classification models to perform hits selection.
        """
        
        # update fname_suffix if applicable
        self.fname_suffix = ''
        if plate is not None:
            self.fname_suffix += f'_{plate}'
        if run is not None:
            self.fname_suffix += f'_{run}'
        if exp_workflow is not None:
            if isinstance(exp_workflow, list):
                self.fname_suffix += f'_{"-".join(exp_workflow)}'
            else:
                self.fname_suffix += f'_{exp_workflow}'
        
        # record Train input
        self.train_input = locals()
        self.train_input.pop('self')
        print('Input to Train method:', self.train_input)
        
        ################################
        # get labeled combined dataset #
        ################################
        
        data = self.GET_COMBINED_DATASETS(
            plate=plate, run=run, exp_workflow=exp_workflow, data_types_to_get=data_types_to_get, model_type_list=self.model_type_list,
              maldi_1_csv=maldi_1_csv, maldi_2_csv=maldi_2_csv, lcms_csv=lcms_csv, lcms_chiral_csv=lcms_chiral_csv,
            get_source_plate_groupings=get_source_plate_groupings, source_plate_groupings=source_plate_groupings, plate_to_grp_idx_mapping=plate_to_grp_idx_mapping,
            split_by_list=split_by_list, get_split_by_list=True, chiral_merge_method='outer',
            get_binary_labels_maldi=True, get_binary_labels_lcms=get_lcms_labels, calculate_conversion_enantioselectivity=True)
        
            
        #####################
        ### train  models ###
        #####################
        
        if not skip_model_training:
            self.models = {}
            self.model_metrics = {}
            for model_type in self.model_type_list:

                ## CLASSIFICATION MODELS
                if model_type.find('classification') > 0:
                    X_colnames = self.maldi_classification_features
                    y_colname = 'true_score_binary'

                    if model_type == 'maldi1sample_classification':
                        data_type = 'maldi_detections_1'
                    elif model_type == 'maldi2sample_classification':
                        data_type = 'maldi_lcms' if 'maldi_lcms' in data else 'maldi_detections_2'

                    for col_suffix in self.split_by_list:
                        model_name = f'{model_type}*{col_suffix}' if col_suffix != '' else model_type
                        print(f'Training {model_name} models with {data_type} data.')
                        # specify whether to use lcms-labeled column or not by appending to col_suffix
                        self.models[model_name], self.model_metrics[model_name] = self.train_maldi2maldi_binary(data[data_type],
                                                                                             X_colnames, y_colname, col_suffix)

                ## REGRESSION MODELS
                elif model_type.find('regression') > 0:

                    data_type = 'maldi_lcms'
                    X_colnames = self.maldi_classification_features
                    y_colname = self.lcms_output_feature

                    # subselect data for regression training
                    if len(self.hit_filter_dict) > 0:
                        df_maldi_lcms_selected = data[data_type][eval(get_metric_filter_str(self.hit_filter_dict, joiner=' | '))]
                        print(f'Filtered regression dataset using hit_filter_dict. Number of samples: {len(df_maldi_lcms_selected)}')
                    else:
                        df_maldi_lcms_selected = data[data_type]

                    for col_suffix in self.split_by_list:
                        model_name = f'{model_type}*{col_suffix}' if col_suffix != '' else model_type
                        print(f'Training {model_name} models with {data_type} data.')
                        self.models[model_name], self.model_metrics[model_name]  = self.train_maldi2lcms_nonbinary(
                            df_maldi_lcms_selected, X_colnames, y_colname, col_suffix)
            
        return data
 
    
    def PREDICT(self, plate=None, run=None, exp_workflow=None, lcms_method='',
                data_types_to_get=['maldi_detections_1', 'maldi_detections_2', 'lcms_detections', 'lcms_detections_chiral'],
                maldi_1_csv=None, maldi_2_csv=None, lcms_csv=None, lcms_chiral_csv=None,
                get_source_plate_groupings=False, source_plate_groupings={}, plate_to_grp_idx_mapping={}, split_by_list_PREDICT=None,
                generate_heatmaps=True, generate_boxplots=True
                ):
        """
        Use trained models to perform hits selection on a plate of data.
        """
        # update fname_suffix if applicable
        self.fname_suffix = ''
        if plate is not None:
            self.fname_suffix += f'_{plate}'
        if run is not None:
            self.fname_suffix += f'_{run}'
        if exp_workflow is not None:
            if isinstance(exp_workflow, list):
                self.fname_suffix += f'_{"-".join(exp_workflow)}'
            else:
                self.fname_suffix += f'_{exp_workflow}'
        
        ##########################################
        # get split_by_list_PREDICT & plate_list #
        ##########################################

        if split_by_list_PREDICT is None:
            self.split_by_list_PREDICT = self.split_by_list
        else:
            self.split_by_list_PREDICT = split_by_list_PREDICT
        
        ################################
        # get labeled combined dataset #
        ################################
        
        data = self.GET_COMBINED_DATASETS(
            plate=plate, run=run, exp_workflow=exp_workflow, data_types_to_get=data_types_to_get, model_type_list=self.model_type_list,
            maldi_1_csv=maldi_1_csv, maldi_2_csv=maldi_2_csv, lcms_csv=lcms_csv, lcms_chiral_csv=lcms_chiral_csv,
            get_source_plate_groupings=get_source_plate_groupings, source_plate_groupings=source_plate_groupings, plate_to_grp_idx_mapping=plate_to_grp_idx_mapping,
            split_by_list=self.split_by_list_PREDICT, get_split_by_list=False, chiral_merge_method='outer',
            get_binary_labels_maldi=True, get_binary_labels_lcms=False, calculate_conversion_enantioselectivity=True)

        ###################################
        ###### predict on maldi data ######
        ###################################
        
        # update predictions on datasets
        for model_type in self.model_type_list:
            
            # classification
            if model_type.find('classification') > 0:
                X_colnames = self.maldi_classification_features
                y_colname = None
                
                if model_type == 'maldi1sample_classification':
                    data_type = 'maldi_detections_1'
                elif model_type == 'maldi2sample_classification':
                    data_type = 'maldi_lcms'
                    
                for col_suffix in self.split_by_list_PREDICT:
                    substrate_barcode = col_suffix
                    model_name = f'{model_type}*{col_suffix}' if col_suffix != '' else model_type
                    models = self.models[model_name]
                    y_pred, y_pred_variance = self.predict_maldi2maldi_binary(models, data[data_type],
                                                                              X_colnames, y_colname, col_suffix)
                    print('Average predicted label:', round(np.nanmean(y_pred),4), \
                          '; Average prediction variance:', round(np.nanmean(y_pred_variance),4))
                    
                    # update values in dataframe
                    col_suffix = col_suffix if col_suffix=='' else f'_{col_suffix}'
                    data[data_type].loc[:,f'predicted_binary_score{col_suffix}'] = list(y_pred)
                    data[data_type].loc[:,f'predicted_binary_score_variance{col_suffix}'] = list(y_pred_variance)
                    print(f'Used {model_name} models to perform prediction on {data_type} data.')
                    
            # regression
            elif model_type.find('regression') > 0:
                data_type = 'maldi_lcms'
                X_colnames = self.maldi_classification_features
                y_colname = None
            
                for col_suffix in self.split_by_list_PREDICT:
                    substrate_barcode = col_suffix
                    model_name = f'{model_type}*{col_suffix}' if col_suffix != '' else model_type
                    models = self.models[model_name]
                    y_pred, y_pred_variance = self.predict_maldi2lcms_nonbinary(models, data[data_type],
                                                                                X_colnames, y_colname, col_suffix)
                    print('Average predicted label:', round(np.nanmean(y_pred),4), \
                          '; Average prediction variance:', round(np.nanmean(y_pred_variance),4))
                    
                    # update values in dataframe
                    col_suffix = col_suffix if col_suffix=='' else f'_{col_suffix}'
                    data[data_type].loc[:,f'predicted_nonbinary_score{col_suffix}'] = list(y_pred)
                    data[data_type].loc[:,f'predicted_nonbinary_score_variance{col_suffix}'] = list(y_pred_variance)
                    print(f'Used {model_name} models to perform prediction on maldi_detections_2 data.')

        
        ########################################
        ###### get enzyme_analytics_df ######
        ########################################
        self.enzyme_analytics_df = self.GET_ENZYME_DATASET(data, split_by_list=self.split_by_list_PREDICT, predict_labels=True)
                    
        ###########################
        ### save datasets in S3 ###
        ###########################
        # predicted datasets
        for data_type, df in data.items():
            s3_df2csv(df, self.s3_bucket, f'{self.s3_subfolder}{data_type}{self.fname_suffix}_PRED.csv')
            print(f'Saved {data_type} dataset to {self.s3_bucket}{self.s3_subfolder} in S3')
            
        # enzyme analytics table
        s3_df2csv(self.enzyme_analytics_df, self.s3_bucket, f'{self.s3_subfolder}enzyme_analytics{self.fname_suffix}_PRED.csv')
        print(f'Saved enzyme analytics dataset to {self.s3_bucket}{self.s3_subfolder} in S3')
        
        ###################################################
        #### update enzyme_analytics table in Postgres ####
        ###################################################
        enzyme_analytics_df = self.update_analytics_table_postgres(enzyme_analytics_df)
    
            
    def HIT_SELECTION(self, enzyme_analytics_df=None, select_by=['predicted_binary_score', 'predicted_nonbinary_score_(r)'], max_num_hits=100, select_up_to_max=False,
                message_prefix='', plate_type_to_select_from='maldi_plate', generate_worklist=False):
        
        #################################
        ###### worklist generation ######
        #################################
        
        # get input variants
        if enzyme_analytics_df is not None:
            self.enzyme_analytics_df = enzyme_analytics_df
            
        # select variants from experimental controls only, then append experimental variants
        ctrl_variants = self.enzyme_analytics_df[self.enzyme_analytics_df.ctrl_type.isin(['pos', self.neg_ctrltype])].copy()
        exp_variants = self.enzyme_analytics_df[self.enzyme_analytics_df.ctrl_type=='exp'].copy()
        self.selected_variants = self.select_variants(exp_variants, select_by=select_by, max_num_hits=max_num_hits,
                                                      select_up_to_max=select_up_to_max, col_suffix_list=self.split_by_list_PREDICT)
        self.selected_variants = pd.concat([self.selected_variants, ctrl_variants])
        print(f'Final # of variants selected (with controls added): {len(self.selected_variants)}')
        
        # generate worklist on a plate by plate basis
        self.worklist_dict = {}
        for col_suffix in self.split_by_list_PREDICT:
            col_suffix_db = self.col_suffix_db[col_suffix]
            if col_suffix_db.find('(r)') > -1:
                plates = list(set(self.selected_variants[f'{plate_type_to_select_from}{col_suffix_db}']))
                for plate in plates:
                    selected_variants_plate = self.selected_variants[self.selected_variants[f'{plate_type_to_select_from}{col_suffix_db}']==plate]
                    print(f'# of variants selected from {plate}: {len(selected_variants_plate)}')
                    if generate_worklist:
                        self.worklist_dict[col_suffix] = self.generate_worklist(selected_variants_plate,
                                                                           track_dilution_plates=self.track_dilution_plates,
                                                                           col_suffix=col_suffix,
                                                                           plate_type_to_select_from=plate_type_to_select_from)

                        # post updates
                        self.post_message_updates(plate, len(selected_variants_plate), message_prefix)

        # update maldi plate status in maldi_detections_intake_queue
        if plate_type_to_select_from == 'maldi_plate':
            table = 'maldi_detections_intake_queue'
        elif plate_type_to_select_from == 'lcms_plate':
            table = 'lcms_detections_intake_queue'
        for plate in self.plate_list:
            DatabaseInterface(table=table).row_update(['analysis_processed'], [True], ['%s'], {'plate': plate})
    
    
    def update_analytics_table_postgres(self, enzyme_analytics_df, enzyme_analytics_table='combi_analytics_table', update_id=True):
        
        # update id column
        if update_id:
            ids = [f'{w}_{p}_{a}' for w, p, a in zip(enzyme_analytics_df['exp_workflow_barcode'].tolist(), enzyme_analytics_df['lcms_plate_(r)'].tolist(), enzyme_analytics_df['lcms_address_(r)'].tolist())]
            enzyme_analytics_df['id'] = ids
        
        # add missing columns to table
        enzyme_analytics_df = DatabaseInterface(table=enzyme_analytics_table).pad_with_missing_columns(enzyme_analytics_df, combi_analytics_metrics)
        # remove rows which don't have any associated source plates
        enzyme_analytics_df = enzyme_analytics_df[~enzyme_analytics_df['source_plate_(+)'].isnull() | ~enzyme_analytics_df['source_plate_(-)'].isnull() | ~enzyme_analytics_df['source_plate_(r)'].isnull()]
        
        # update enzyme_analytics_table records in S3 by unique source address, source plate combinations
        if self.update_enzyme_analytics_table:
            source_plate_list = list(set(enzyme_analytics_df['source_plate_(r)'].tolist()))
            print('source_plate_(r) list:', source_plate_list)
            for source_plate in source_plate_list:
                # get data relevant to source plate
                enzyme_analytics_df_plate = enzyme_analytics_df[enzyme_analytics_df['source_plate_(r)']==source_plate].copy()
                # clear previous data for source plate
                db_detections = DatabaseInterface(table=enzyme_analytics_table)
                db_detections.clear_matching_records({'source_plate_(r)': source_plate})
                # update with new data
                db_detections.table_add(enzyme_analytics_df_plate)
                print(f'Updated postgres enzyme_analytics_table for all (+), (-), (r) samples associated with the variants on source plate {source_plate}.')
        return enzyme_analytics_df
