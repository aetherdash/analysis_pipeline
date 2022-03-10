import os, sys
sys.path.append('..')
import random
import pandas as pd
import numpy as np
import pickle
import json
import itertools
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

from analytics_utils.database_access.table_properties import *
from analytics_utils.database_access.db_interface import DatabaseInterface
from analytics_utils.database_access.s3_interface import download_from_s3, upload_to_s3, s3_imgupload, s3_df2csv
from analytics_utils.analysis_tools.analysis_utils import get_ctrl_vals, get_lcms_qc_derivedmetrics_checks, calculate_conversion_enantioselectivity_C18, calculate_stats_for_array, calculate_stats_for_dataframe, get_unique_variants_in_dataframe
from analytics_utils.visualization_tools.visualization_utils import generate_plots, df2array_dict, plot_boxplot
from analytics_utils.lims_tools.lims_utils import lims_post_matplotlib_img
        
class LcmsChiralAnalytics:
    
    def __init__(self,
                 s3_bucket='ml-models-registry',
                 s3_subfolder='exerevnetes-preprocessing-models/ProjectX_UnitX/',
                 fname_suffix='',
                 posctrl_enzyme_barcode = None,
                 maldi_label_condition = {},
                 hit_filter_dict = {},
                 plot_enantiomer_lcms_metrics = None,
                 col_suffix_db = {'':'_(r)', 'CMP60354':'_(r)', 'CMP60403':'_(-)', 'CMP60404':'_(+)'},
                 neg_ctrltype='EV',
                 local_dir='../DATASETS/',
                 update_plate_analytics_table=True,
                 update_enzyme_analytics_table=True,
                 track_dilution_plates=False,
                 save_plots_to_lims=True,
                 get_2sample_from_1sample_dataset=True
                ):

        self.s3_bucket = s3_bucket
        self.s3_subfolder = s3_subfolder
        self.fname_suffix = fname_suffix
        self.local_dir = local_dir
        self.metricsetname = 'LcmsChiral'
        self.cols_to_get = {'lcms_detections_chiral': lcms_detections_chiral_columns_to_get}
        self.neg_ctrltype = neg_ctrltype
        self.posctrl_enzyme_barcode = posctrl_enzyme_barcode
        self.plot_enantiomer_lcms_metrics = plot_enantiomer_lcms_metrics
        self.col_suffix_db = col_suffix_db
        self.save_plots_to_lims = save_plots_to_lims
        self.get_2sample_from_1sample_dataset = get_2sample_from_1sample_dataset
        self.ctrl_to_color_mapping = {'exp':'r', 'pos':'g', 'EV':'k'}
        self.stats = ['mean','std','median','iqr','iqr_norm','cv', 'max', 'min', 'range', 'range_norm']
    
    
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
            
        return df
    
    
    def get_derived_metrics(self, df):
        
        df_wderived_metrics = df.copy()
        conversion, ee_plusoverminus, eeratio_plusoverminus = calculate_conversion_enantioselectivity_chiral(df_wderived_metrics, ee_calc_method=3)
        for s in conversion:
            df_wderived_metrics.loc[:,f'measured_conversion{s}_chiral'] = conversion[s]
        df_wderived_metrics.loc[:,'measured_enantiomeric_excess_(+over-)_chiral'] = ee_plusoverminus
        df_wderived_metrics.loc[:,'measured_enantiomeric_ratio_(+over-)_chiral'] = eeratio_plusoverminus

        return df_wderived_metrics
    
    
    def get_qc_metrics(self):
        return None

     
    def get_group_analytics(self, df, groupby, grp_list, metricset, neg_ctrltype='EV', metricsetname=None,
                           filter_posctrl=('ctrl_type', ['pos']), filter_negctrl=('ctrl_type', ['EV']), filter_expctrl=('ctrl_type', ['exp']),
                            calc_FIOP_bygrp=False, get_ctrls_from_vals=False):

        if metricsetname is not None:
            self.metricsetname = metricsetname

        ## calculate metrics by metricset ##
        cols_to_get = ['ctrl_type', 'enzyme_barcode'] + [g for g in groupby if g not in ['ctrl_type', 'enzyme_barcode']] + metricset
        df.loc[:, metricset] = df[metricset].fillna(value=np.nan)
        df_exp = df.loc[df[filter_expctrl[0]].isin(filter_expctrl[1]), cols_to_get].copy()
        df_pos = df.loc[df[filter_posctrl[0]].isin(filter_posctrl[1]), cols_to_get].copy()
        df_neg = df.loc[df[filter_negctrl[0]].isin(filter_negctrl[1]), cols_to_get].copy()
        df_ctrl = pd.concat([df_pos, df_neg])
        df_ctrls_dict = {'pos':pd.DataFrame(columns=cols_to_get), self.neg_ctrltype: pd.DataFrame(columns=cols_to_get)}
 
        # initialize containers
        stats_bymetricset_CTRL_GRP = {}

        # get GROUP stats for POS, NEG ctrls
        for grp in grp_list:
            ## get data for POS + NEG ctrls for group
            df_ctrl_grp_dict = {}
            if get_ctrls_from_vals:
                df_ctrl_grp = df_ctrl.loc[(df_ctrl[groupby[0]]==grp)].copy()
                df_ctrl_grp_dict['pos'], df_ctrl_grp_dict[self.neg_ctrltype] = get_ctrl_vals(df_ctrl_grp, cols_to_get=cols_to_get, col_to_sort_by=metricset[0], get_ctrls_from_vals=True, compare_byctrlval=True, n_max=None, n_min=None, frac_max=0.25, frac_min=0.25, neg_ctrltype=self.neg_ctrltype)
            else:
                df_ctrl_grp_dict['pos'] = df_pos.loc[(df_pos[groupby[0]]==grp)].copy()
                df_ctrl_grp_dict[self.neg_ctrltype] = df_neg.loc[(df_neg[groupby[0]]==grp)].copy()

            # update overall stats dict
            df_ctrls_dict['pos'] = df_ctrls_dict['pos'].append(df_ctrl_grp_dict['pos'], ignore_index=True)
            df_ctrls_dict[self.neg_ctrltype] = df_ctrls_dict[self.neg_ctrltype].append(df_ctrl_grp_dict[self.neg_ctrltype], ignore_index=True)
#             print(grp, '# of pos ctrl samples:', len(df_ctrl_grp_dict['pos']), '; # of EV ctrl samples:', len(df_ctrl_grp_dict[self.neg_ctrltype]))
                
            # get GROUP metric stats POS + NEG ctrls
            CTRL_stats_bygrp = {}
            for ctrl_type in ['pos', self.neg_ctrltype]:
                for metric in metricset:
                    df_CTRL_grp = df_ctrl_grp_dict[ctrl_type]
                    CTRL_stats_bygrp = calculate_stats_for_dataframe(df_CTRL_grp, metric, stats_to_get=self.stats, metrics_dict=CTRL_stats_bygrp,
                                                                     metric_prefix=f'{metric}_', metric_suffix=f'_{ctrl_type.upper()}')
            stats_bymetricset_CTRL_GRP.update({grp: CTRL_stats_bygrp})

        # get OVERALL stats for POS, NEG ctrls
        stats_table_bymetricset_CTRL, stats_bymetricset_CTRL_ALL, CTRL_samples_bymetricset = self.get_ctrl_stats_bymetricset(df_ctrls_dict, groupby, metricset)

        # get GROUP metric stats for EXP ctrls
        stats_table_bymetricset_EXP_GRP = []
        for grp in grp_list:

            ## get data for EXP ctrls for grp
            stats_bygrp = {}
            df_grp = df_exp.loc[df_exp[groupby[0]]==grp]
            
            if len(df_grp) > 0:
                # get additional metadata
                for g in groupby:
                    stats_bygrp.update({f'{g}': df_grp.iloc[0][g]})
                # update number of unique variants
                num_unique_variants, _, _ = get_unique_variants_in_dataframe(df_grp, colname='enzyme_barcode', metric_prefix='', metric_suffix='')
                stats_bygrp.update({f'n_var': num_unique_variants})

                for i, metric in enumerate(metricset):
                    if i==len(metricset)-1:
                        n = len(df_grp[~df_grp[metric].isnull()])
                        stats_bygrp.update({f'n_{self.metricsetname}': n})
                    # get pos ctrl metrics for calculating FIOP
                    if calc_FIOP_bygrp:
                        pos_median = stats_bymetricset_CTRL_GRP[grp][f'{metric}_median_POS']
                    else:
                        pos_median = stats_bymetricset_CTRL_ALL['pos'][f'{metric}_median']

                    # get exp stats
                    stats_bygrp = calculate_stats_for_dataframe(df_grp, metric, stats_to_get=self.stats, metric_prefix=f'{metric}_', metrics_dict=stats_bygrp)
                    # get FIOP & update metrics
                    median = stats_bygrp[f'{metric}_median']
                    maximum = stats_bygrp[f'{metric}_max']
                    median_FIOP = float(median/pos_median)
                    max_FIOP = float(maximum/pos_median)
                    stats_bygrp.update({f'{metric}_max_FIOP': max_FIOP})
                    stats_bygrp.update({f'{metric}_median_FIOP': median_FIOP})

                # update exp_stats with POS stats
                stats_bygrp.update(stats_bymetricset_CTRL_GRP[grp])
                stats_table_bymetricset_EXP_GRP.append(stats_bygrp)
                
        stats_table_bymetricset_EXP_GRP = pd.DataFrame(stats_table_bymetricset_EXP_GRP)
                
        # get overall EXP stats
        stats_table_bymetricset_EXP_ALL = self.get_exp_stats_bymetricset(stats_table_bymetricset_EXP_GRP, groupby)
     
        # compile overall stats_table
        stats_table_bymetricset = pd.concat([stats_table_bymetricset_CTRL, stats_table_bymetricset_EXP_ALL, stats_table_bymetricset_EXP_GRP], axis=0)
                
        return stats_table_bymetricset, stats_bymetricset_CTRL_GRP, stats_bymetricset_CTRL_ALL, CTRL_samples_bymetricset
    
    
    def get_ctrl_stats_bymetricset(self, df_ctrls_dict, groupby_list, metricset):
        """
        Get OVERALL (NOT by grp) stats for POS & NEG ctrls
        """
        CTRL_samples_bymetricset = {}
        ctrltype_list = [self.neg_ctrltype, 'pos']
        index_list = [f'all {ctrl_type} ctrls' for ctrl_type in ctrltype_list]
        stats_table_bymetricset_CTRL = []
        stats_bymetricset_CTRL_ALL = {ctrltype:{} for ctrltype in ctrltype_list}
        for ctrl_type, idx in zip(ctrltype_list, index_list):
            df_ctrl = df_ctrls_dict[ctrl_type]
            # get additional metadata
            for g in groupby_list:
                stats_bymetricset_CTRL_ALL[ctrl_type].update({f'{g}': idx})
            for i, metric in enumerate(metricset):
                if i==len(metricset)-1:
                    n = len(df_ctrl[~df_ctrl[metric].isnull()])
                    stats_bymetricset_CTRL_ALL[ctrl_type].update({f'n_{self.metricsetname}': n})
                stats_bymetricset_CTRL_ALL[ctrl_type] = calculate_stats_for_dataframe(df_ctrl, metric, stats_to_get=self.stats, metric_prefix=f'{metric}_', metrics_dict=stats_bymetricset_CTRL_ALL[ctrl_type])
            stats_table_bymetricset_CTRL.append(stats_bymetricset_CTRL_ALL[ctrl_type])
            
            # collect all ctrl samples for particular metric
            CTRL_samples_bymetricset[ctrl_type] = df_ctrl # {self.metricsetname:df_ctrl}
            
        stats_table_bymetricset_CTRL = pd.DataFrame(stats_table_bymetricset_CTRL, index=index_list)
        return stats_table_bymetricset_CTRL, stats_bymetricset_CTRL_ALL, CTRL_samples_bymetricset
    
    def get_exp_stats_bymetricset(self, stats_table_bymetricset_EXP_GRP, groupby_list):
        """
        Get average of grps stats for EXP ctrls
        """
        # get numeric cols
        stats_table_bymetricset_EXP_GRP_numeric = stats_table_bymetricset_EXP_GRP.select_dtypes(include=['int','float'])
        numeric_cols = stats_table_bymetricset_EXP_GRP_numeric.columns.tolist()
        
        # compute average stats by group 
        stats_table_bymetricset_EXP_ALL = pd.DataFrame(stats_table_bymetricset_EXP_GRP[numeric_cols].mean(axis=0, skipna=True)).transpose()
        stats_table_bymetricset_EXP_ALL.index = [f'all exp ctrls (by {groupby_list[0]})']
        return stats_table_bymetricset_EXP_ALL
    
    
    def get_overall_visualizations(self, df, columns_to_plot, groupby):
        fig_list = plot_boxplot(df, columns_to_plot, groupby=[groupby], fname_prefix='', fname_suffix=f'_{groupby}', plate_barcode=None, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder, save_plot_to_s3=True, save_plot_to_lims=False)
        return fig_list
    
    def get_plate_visualizations(self, df, columns_to_plot, groupby):
        fig_list = plot_boxplot(df, columns_to_plot, groupby=[groupby, 'ctrl_type'], fname_prefix='', fname_suffix=f'_{groupby}', plate_barcode=None, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder, save_plot_to_s3=True, save_plot_to_lims=False)
        return fig_list
    
    def get_library_visualizations(self, df, stats_table_bymetricset, columns_to_plot, groupby='library_barcode', top_n_libraries=100, sort_by=('measured_enantiomeric_excess_(+over-)', 'median')):
        from collections import OrderedDict

        # get top n libraries
        metric, stat = sort_by[0], sort_by[1]
        stats_table_filt = stats_table_bymetricset.iloc[3:].sort_values(f'{metric}_{stat}', ascending=False)
        if top_n_libraries is not None:
            stats_table_filt = stats_table_filt.iloc[:min(top_n_libraries, len(stats_table_filt))]
        lib_list = stats_table_filt[groupby].tolist()
        
        # get groups based on top libraries
        df_orderedgrps_dict = {column: OrderedDict() for column in columns_to_plot}
        for lib in lib_list: 
            for column in columns_to_plot:
                df_orderedgrps_dict[column][lib] = df.loc[df[groupby]==lib, column]

        # get boxplots
        fig_list = []
        for column in columns_to_plot:
            figs = plot_boxplot(pd.DataFrame(df_orderedgrps_dict[column]), [column], groupby=None, fname_prefix='', fname_suffix=f'_library', plate_barcode=None, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder, save_plot_to_s3=True, save_plot_to_lims=False)
            fig_list += figs
        return fig_list
    
    def get_unit_visualizations(self, df, columns_to_plot, groupby):
        fig_list = plot_boxplot(df, columns_to_plot, groupby=[groupby, 'ctrl_type'], fname_prefix='', fname_suffix=f'_{groupby}', plate_barcode=None, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder, save_plot_to_s3=True, save_plot_to_lims=False)
        return fig_list
    
#     def plot_scatter_chiral(self, df, groupby=None, grp_list=None):
#         """
#         Plot chiral enantiomer activity or conversion as scatter plot
#         """
#        # CONVERSION + VS -
#         plot_lcms_metrics_for_ctrls(df,
#                                     x_suffix=f'(-)',
#                                     y_suffix=f'(+)',
#                                     ctrls_to_plot=['exp','pos',self.neg_ctrltype],
#                                     metric_list=['measured_conversion'])
#         # upload plot to S3
#         s3_imgupload(self.s3_bucket, f"{self.s3_subfolder}lcms_enantiomer_metrics_CONVERSION{self.fname_suffix}.jpg")
#         plt.show()
