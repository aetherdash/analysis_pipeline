import os, sys
sys.path.append('..')
import pandas as pd
import numpy as np
import panel as pn
import pickle
from collections import OrderedDict
from datetime import datetime
import matplotlib.pyplot as plt

from analytics_utils.database_access.table_properties import *
from analytics_utils.database_access.sql_utils import get_metric_filter_str
from analytics_utils.database_access.db_interface import DatabaseInterface
from analytics_utils.database_access.s3_interface import download_from_s3, upload_to_s3, s3_imgupload, s3_imgdownload, s3_df2csv, s3_csv2df
from analytics_utils.maldi_processing.ttest import get_peak_vals, get_ttests
from analytics_utils.analysis_tools.analysis_utils import get_ctrl_vals, get_lcms_qc_derivedmetrics_checks, calculate_conversion_enantioselectivity_C18, calculate_conversion_enantioselectivity_chiral, calculate_stats_for_array, calculate_stats_for_dataframe, get_unique_variants_in_dataframe
from analytics_utils.visualization_tools.visualization_utils import generate_plots, df2array_dict, plot_scatter_combined_metrics, plot_histogram_combined_metrics, plot_boxplot
from analytics_utils.lims_tools.lims_utils import lims_post_matplotlib_img

# Analytics pipeline classes
from .lcms_C18_analytics import LcmsC18Analytics
from .lcms_chiral_analytics import LcmsChiralAnalytics
        
class AnalysisPipeline:
    
    def __init__(self,
                 s3_bucket='ml-analytics-file-store',
                 s3_subfolder='ProjectX/',
                 fname_suffix='',
                 posctrl_enzyme_barcode = None,
                 col_suffix_db = {'':'_(r)', 'CMP60354':'_(r)', 'CMP60403':'_(-)', 'CMP60404':'_(+)'},
                 neg_ctrltype='EV',
                 local_dir='../DATASETS/',
                 update_plate_analytics_table=True,
                 track_dilution_plates=False,
                 save_plots_to_s3=True,
                 save_plots_to_lims=False,
                 sort_by=None,
                 metric_dict={'LcmsC18':['pellet_OD', 'measured_conversion_(r)'], 'LcmsChiral':['measured_enantiomeric_excess_(+over-)']},
                 get_dashboard_panel=False,
                 hit_filter_dict=None
                ):

        self.s3_bucket = s3_bucket
        self.s3_subfolder = s3_subfolder
        self.fname_suffix = fname_suffix
        self.local_dir = local_dir
        self.cols_to_get = {
            'maldi_detections_1': maldi_detections_columns_to_get,
            'maldi_detections_2': maldi_detections_columns_to_get,
            'lcms_detections': lcms_detections_columns_to_get,
            'lcms_detections_chiral': lcms_detections_chiral_columns_to_get,
            'combi_analytics_table': combi_analytics_metrics
                           }
        self.posctrl_enzyme_barcode = posctrl_enzyme_barcode
        self.neg_ctrltype = neg_ctrltype
        self.save_plots_to_s3 = save_plots_to_s3
        self.save_plots_to_lims = save_plots_to_lims
        self.analytics_classes = {
            'LcmsC18':LcmsC18Analytics(neg_ctrltype=self.neg_ctrltype, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder, hit_filter_dict=hit_filter_dict),
            'LcmsChiral':LcmsChiralAnalytics(neg_ctrltype=self.neg_ctrltype, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder)
        }
        self.metric_dict = metric_dict
        self.metric_list = sum(self.metric_dict.values(), [])
        self.scatterplot_axes = {
            'prod(+)C18_vs_prod(-)C18':(['measured_nonbinary_score_(-)','measured_nonbinary_score_(+)'], ['prod_conc_(-)_C18', 'prod_conc_(+)_C18']),
            'EE_vs_Conversion':(['measured_conversion_(r)','measured_enantiomeric_excess_(+over-)'], ['RacemicConversion', 'ChiralEE']),
        }
        # self.stats = ['median','mean','max', 'min','cv','std','iqr','iqr_norm','range', 'range_norm']
        # self.suffix_list = [f'_{stat}_{ctrl_type}' for ctrl_type in [self.neg_ctrltype, 'POS'] for stat in ['median','mean','cv','std','iqr']] + \
        # ['_max_FIOP', '_median_FIOP'] + [f'_{stat}' for stat in self.stats+['nonhitrate']]
        self.stats = ['median','mean','max', 'cv','std',]
        self.suffix_list = [f'_{stat}_{ctrl_type}' for ctrl_type in [self.neg_ctrltype, 'POS'] for stat in ['median','mean','cv','std']] + \
        ['_max_FIOP', '_median_FIOP'] + [f'_{stat}' for stat in self.stats+['nonhitrate']]
        self.suffix_list_to_display = ['_max_FIOP', '_median_FIOP'] + [f'_{stat}' for stat in ['median','mean','cv','nonhitrate']]
        self.sort_by = sort_by
        self.get_dashboard_panel = get_dashboard_panel
        self.dashboard_contents = []
    
        
    
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
    
    
    def get_derived_metrics(self, df, derived_metrics_to_get=[]):
        
        for metricsetname, metricset in self.metric_dict.items():
            if metricsetname in derived_metrics_to_get:
                df = self.analytics_classes[metricsetname].get_derived_metrics(df)
                
        self.df = df
        return df
    
    
    def get_qc_metrics(self, df, derived_metrics_to_get=[]):
        
        for metricsetname, metricset in self.metric_dict.items():
            if metricsetname in derived_metrics_to_get:
                df = self.analytics_classes[metricsetname].get_derived_metrics(df)
                
        return df
        
    
    def standardize_dataset_columns(self, df, data_type='combi_analytics_table', update_analytics_table=False):
                      
        if data_type == 'lcms_detections':
            cols_to_get = [
                'exp_workflow_barcode',
                'exp_workflow_name',
                'proj_barcode',
                'proj_name',
                'run',
                'plate',
                'address',
                'source_plate',
                'source_address',
                'ctrl_type',
                'exp_condition',
                'enzyme_barcode',
                'sequence',
                'mutations',
                'hamming',
                'reference_enzyme',
                'substrate_barcode',
                'substrate_smiles',
                'substrate_concentration',
                'prod_conc_lcms_actual',
                'measured_conversion_(r)',
                'measured_enantiomeric_excess_(+over-)',
                'sum_conc_lcms_actual'] + \
            ['seed_address', 'seed_address_alphanumeric', 'seed_plate', 'library_barcode', 'library_ref', 'library_description', 'seed_plate_time'] + \
            pelletqc_metadata
        
            cols_to_rename = {
                'run': 'lcms_C18_run',
                'address':'lcms_address_(r)',
                'plate':'lcms_plate_(r)',
                'source_plate':'source_plate_(r)',
                'source_address':'source_address_(r)',
                'substrate_barcode':'substrate_barcode_(r)',
                'substrate_smiles':'substrate_smiles_(r)',
                'substrate_concentration': 'substrate_concentration_(r)',
                'prod_conc_lcms_actual':'measured_nonbinary_score_(r)',
                'sum_conc_lcms_actual': 'measured_nonbinary_sum_(r)'
            }
            
        elif data_type=='lcms_detections_all':
            cols_to_get = [
                'exp_workflow_barcode', 
                'exp_workflow_name', 
                'proj_barcode', 
                'proj_name', 
                'run_CMP60354', 
                'run', 
                'exp_condition',
                'plate_CMP60354', 
                'address', 
                'source_plate_CMP60354', 
                'source_address_CMP60354', 
                'ctrl_type', 
                'enzyme_barcode', 
                'sequence', 
                'mutations', 
                'hamming', 
                'reference_enzyme', 
                'substrate_barcode_CMP60354', 
                'substrate_smiles_CMP60354', 
                'substrate_concentration_CMP60354', 
                'prod_conc_lcms_actual_CMP60354',
                'sum_conc_lcms_actual_CMP60354',
                'sub_conc_lcms_actual_(-)', 
                'sub_conc_lcms_actual_(+)', 
                'measured_binary_score_(r)_lcms', 
                'measured_conversion_(r)_C18', 
                'measured_enantiomeric_excess_(+over-)_chiral', 
                'measured_enantiomeric_ratio_(+over-)_chiral'] + \
            ['seed_address', 'seed_address_alphanumeric', 'seed_plate', 'library_barcode', 'library_ref', 'library_description', 'seed_plate_time'] + \
            pelletqc_metadata
            
            cols_to_rename = {
                'run_CMP60354': 'lcms_C18_run', 
                'run':'lcms_chiral_run', 
                'address':'lcms_address_(r)', 
                'plate_CMP60354':'lcms_plate_(r)', 
                'source_plate_CMP60354':'source_plate_(r)', 
                'source_address_CMP60354':'source_address_(r)', 
                'substrate_barcode_CMP60354':'substrate_barcode_(r)', 
                'substrate_smiles_CMP60354':'substrate_smiles_(r)', 
                'substrate_concentration_CMP60354': 'substrate_concentration_(r)', 
                'prod_conc_lcms_actual_CMP60354':'measured_nonbinary_score_(r)', 
                'sum_conc_lcms_actual_CMP60354':'measured_nonbinary_sum_(r)', 
                'sub_conc_lcms_actual_(-)': 'measured_nonbinary_score_(-)', 
                'sub_conc_lcms_actual_(+)': 'measured_nonbinary_score_(+)', 
                'measured_binary_score_(r)_lcms':'measured_binary_score_(r)', 
                'measured_conversion_(r)_C18':'measured_conversion_(r)', 
                'measured_enantiomeric_excess_(+over-)_chiral':'measured_enantiomeric_excess_(+over-)', 
                'measured_enantiomeric_ratio_(+over-)_chiral':'measured_enantiomeric_ratio_(+over-)'
            } 
        
        df = df[[c for c in cols_to_get if c in df]].copy().rename(columns=cols_to_rename)
        if update_analytics_table:
            df = self.update_analytics_table_postgres(df, update_analytics_table=update_analytics_table)
        self.df = df
        return df
        
        
    def GET_OVERALL_ANALYTICS(self, df, groupby=['ctrl_type'], display_table=True, plot_scatterplot=True, plot_histogram=True, calc_FIOP_bygrp=False, get_ctrls_from_vals=False):
        print('******************************')
        print('Computing Overall Analytics...')
        print('******************************')

        grp_list = sorted([l for l in list(set(df[groupby[0]])) if isinstance(l, str)])
        stats_table_dict, CTRL_stats, metricset_plots = {}, {}, {}
        CTRL_samples = {'pos':{}, self.neg_ctrltype:{}}

        for metric in self.metric_list:
            df.loc[:, metric] = df[metric].fillna(value=np.nan)

        ########################################################################
        # Get Analytic Metrics by MetricSet (e.g. Maldi, LcmsC18, Lcms Chiral) #
        ########################################################################
        
        for metricsetname, metricset in self.metric_dict.items():
            print(metricsetname, metricset)
            metricset_analysis = self.analytics_classes[metricsetname]
            stats_table_bymetricset, stats_bymetricset_CTRL_GRP, stats_bymetricset_CTRL_ALL, CTRL_samples_bymetricset = metricset_analysis.get_group_analytics(df, groupby, grp_list, metricset, calc_FIOP_bygrp=calc_FIOP_bygrp, get_ctrls_from_vals=get_ctrls_from_vals)
            
            # update joint containers
            stats_table_dict.update({metricsetname:stats_table_bymetricset})
            CTRL_stats.update(stats_bymetricset_CTRL_GRP)
            CTRL_samples['pos'].update(CTRL_samples_bymetricset['pos'])
            CTRL_samples[self.neg_ctrltype].update(CTRL_samples_bymetricset[self.neg_ctrltype])
            
            ## get plots by metricset ##
            metricset_plots[metricsetname] = metricset_analysis.get_overall_visualizations(df, metricset, groupby[0])
            
        #################################
        ## get plots across metricsets ##
        #################################
        metricset, metricnames = self.scatterplot_axes['EE_vs_Conversion']
        scatterplot, histogram = None, None
        # scatter plot
        if plot_scatterplot:
            # Chiral EE vs. Racemic Converson
            metricset, metricnames = self.scatterplot_axes['EE_vs_Conversion']
            scatterplot = plot_scatter_combined_metrics(df, metricset=metricset, metricnames=metricnames, save_plot_to_s3=True, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder, fname_suffix=self.fname_suffix)
            
        # histogram
        if plot_histogram:
            histogram = plot_histogram_combined_metrics(df, CTRL_samples, metricset, metricnames, groupby[0], grp_list, groupbyname='ctrl_type', save_plot_to_s3=True, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder, fname_suffix=self.fname_suffix)

        ################################
        ## Compile Metric Stats Table ##
        ################################
        cols_to_get = groupby + ['n_var'] + [f'n_{metricset}' for metricset in self.metric_dict] + [f'{metric}{suffix}' for metric in self.metric_list for suffix in self.suffix_list]
        cols_to_display = groupby + ['n_var'] + [f'n_{metricset}' for metricset in self.metric_dict] + [f'{metric}{suffix}' for metric in self.metric_list for suffix in self.suffix_list_to_display]
        stats_table = self.get_analytics_table(stats_table_dict, cols_to_get, groupbyname='ctrl_type', sort_by=self.sort_by, cols_to_display=cols_to_display)
        
        # Total number of substrates (GS)
        # Graphical representation of Activity and EE across time (CUMULATIVE) overall [PJ]
        
        # #########################
        # ## Create Panel Column ##
        # #########################
        analytics_panel = []
        if self.get_dashboard_panel:
            text = 'OVERALL ANALYTICS'
            analytics_panel = self.create_panel_columns(text=text, metricset_plots=metricset_plots, scatterplot=scatterplot, histogram=histogram, stats_table=stats_table, groupby='overall', cols_to_display=cols_to_display)
        
        print('Obtained overall analytics metrics. \n')
        return stats_table, analytics_panel 

        
    def GET_VARIANT_ANALYTICS(self, df, groupby=['enzyme_barcode', 'mutations', 'library_ref', 'substrate_concentration_(r)', 'exp_workflow_barcode', 'lcms_plate_(r)', 'lcms_address_(r)'], sort_by=('measured_conversion_(r)', 'median'), display_table=True, plot_scatterplot=False, plot_histogram=True, calc_FIOP_bygrp=False, get_ctrls_from_vals=False, table_suffix=''):
        
        print('******************************')
        print('Computing Variant Analytics...')
        print('******************************')
        
        grp_list = sorted([l for l in list(set(df[groupby[0]])) if isinstance(l, str)])
        stats_table_dict, CTRL_stats = {}, {}
        CTRL_samples = {'pos':{}, self.neg_ctrltype:{}}
        for metric in self.metric_list:
            df.loc[:, metric] = df[metric].fillna(value=np.nan)

        ########################################################################
        # Get Analytic Metrics by MetricSet (e.g. Maldi, LcmsC18, Lcms Chiral) #
        ########################################################################
        
        for metricsetname, metricset in self.metric_dict.items():
            print(metricsetname, metricset)
            metricset_analysis = self.analytics_classes[metricsetname]
            stats_table_bymetricset, stats_bymetricset_CTRL_GRP, stats_bymetricset_CTRL_ALL, CTRL_samples_bymetricset = metricset_analysis.get_group_analytics(df, groupby, grp_list, metricset, calc_FIOP_bygrp=calc_FIOP_bygrp, get_ctrls_from_vals=get_ctrls_from_vals)
            
            # update joint containers
            stats_table_dict.update({metricsetname:stats_table_bymetricset})
            CTRL_stats.update(stats_bymetricset_CTRL_GRP)
            CTRL_samples['pos'].update(CTRL_samples_bymetricset['pos'])
            CTRL_samples[self.neg_ctrltype].update(CTRL_samples_bymetricset[self.neg_ctrltype])
        
        ################################
        ## Compile Metric Stats Table ##
        ################################
        cols_to_get = groupby + [f'n_{metricset}' for metricset in self.metric_dict] + [f'{metric}{suffix}' for metric in self.metric_list for suffix in self.suffix_list]
        cols_to_display = groupby + [f'n_{metricset}' for metricset in self.metric_dict] + [f'{metric}{suffix}' for metric in self.metric_list for suffix in self.suffix_list_to_display]
        stats_table = self.get_analytics_table(stats_table_dict, cols_to_get, groupbyname='enzyme', sort_by=sort_by, table_suffix=table_suffix, cols_to_display=cols_to_display)
            
        #################################
        ## get plots across metricsets ##
        #################################
        metricset, metricnames = self.scatterplot_axes['EE_vs_Conversion']
        scatterplot, histogram = None, None
        # scatter plot
        if plot_scatterplot:
            scatterplot = plot_scatter_combined_metrics(df, metricset, metricnames, groupby[0], grp_list, groupbyname='library', save_plot_to_s3=True, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder, fname_suffix=self.fname_suffix)
        # histogram
        if plot_histogram:
            histogram = plot_histogram_combined_metrics(stats_table.iloc[3:], CTRL_samples, metricset, metricnames, groupby=groupby[0], groupbyname='enzyme', plot_exp_metric_stats=['median', 'iqr'], showplot=False, neg_ctrltype='EV', fig=None, ax=None, save_plot_to_s3=True, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder, fname_suffix=self.fname_suffix)
        
        # #########################
        # ## Create Panel Column ##
        # #########################
        analytics_panel = []
        if self.get_dashboard_panel:
            text = 'VARIANT ANALYTICS'
            analytics_panel = self.create_panel_columns(text=text, scatterplot=scatterplot, histogram=histogram, stats_table=stats_table, groupby='variant', cols_to_display=cols_to_display)
            
        print('Obtained variant analytics metrics. \n')
        return stats_table, analytics_panel
    
    
    def get_top_variants(self, stats_table, df=None, thres_dict={'measured_conversion_(r)':(0.4,1,'median'), 'measured_enantiomeric_excess_(+over-)':(0.4,1,'median')}, top_n_variants=None, sort_by=('measured_conversion_(r)', 'median'), prioritize_filter='n', get_replicates_data=False, display_table=True, table_suffix='', groupby=['enzyme_barcode', 'mutations', 'library_ref', 'substrate_concentration_(r)', 'exp_workflow_barcode', 'lcms_plate_(r)', 'lcms_address_(r)']):
        
        print('**********************************')
        print('Computing Top Variant Analytics...')
        print('**********************************')
        
        # only data with valid enzyme barcode 'ENZxxxxx'
        enz_list = stats_table.enzyme_barcode.tolist()
        idx_enz = [i for i, enz in enumerate(enz_list) if isinstance(enz, str) and enz.upper()[:3]=='ENZ']
        stats_table_filt = stats_table.iloc[idx_enz].copy()
            
        # get top n variants
        if top_n_variants is not None:
            stats_table_filt = stats_table_filt.iloc[:min(top_n_variants, len(stats_table_filt))]
        
        # get data filtered by threshold conditions
        if prioritize_filter=='metric_thres':
            for metric, (thres, inequality, stat) in thres_dict.items():
                if inequality == 1: # less than
                    stats_table_filt = stats_table_filt.loc[stats_table_filt[f'{metric}_{stat}']>thres]
                elif inequality == 0:
                    stats_table_filt = stats_table_filt.loc[stats_table_filt[f'{metric}_{stat}']==thres]
                elif inequality == -1:
                    stats_table_filt = stats_table_filt.loc[stats_table_filt[f'{metric}_{stat}']<thres]

        # get dataframe with replicates data
        if get_replicates_data:
            top_n_variants = stats_table_filt.enzyme_barcode.tolist()
            df_reps_topvariants = pd.DataFrame(columns=df.columns.tolist())
            
            # iterate through selected variants and append all reps seen in order
            for var in top_n_variants:
                df_reps_topvariants = df_reps_topvariants.append(df.loc[df.enzyme_barcode==var].sort_values(sort_by[0], ascending=False))
        else:
            df_reps_topvariants = None
            
        # reset index
        stats_table_filt = stats_table_filt.reset_index(drop=True)
        df_reps_topvariants = df_reps_topvariants.reset_index(drop=True)
        
        # display tables
        if display_table and not self.get_dashboard_panel:
            print(f'Top variant stats (n={len(stats_table_filt)})')
            display(stats_table_filt)
            print('Top variant detections data')
            display(df_reps_topvariants)

        # save tables to S3
        s3_df2csv(stats_table_filt, self.s3_bucket, f'{self.s3_subfolder}top_enzyme_METRICS{self.fname_suffix}{table_suffix}.csv')
        s3_df2csv(df_reps_topvariants, self.s3_bucket, f'{self.s3_subfolder}top_enzyme_DATA{self.fname_suffix}{table_suffix}.csv')
        
        # #########################
        # ## Create Panel Column ##
        # #########################
        cols_to_display = groupby + [f'n_{metricset}' for metricset in self.metric_dict] + [f'{metric}{suffix}' for metric in self.metric_list for suffix in self.suffix_list_to_display]
        
        analytics_panel = []
        if self.get_dashboard_panel:
            text = 'TOP VARIANT ANALYTICS'
            analytics_panel = self.create_panel_columns(text=text, stats_table=None, data_reps=df_reps_topvariants, groupby='top_enzyme', cols_to_display=cols_to_display)
        
        print('Obtained top variant analytics metrics. \n')
        return stats_table_filt, df_reps_topvariants, analytics_panel
        
        
    def GET_PLATE_ANALYTICS(self, df, display_table=True, groupby=['lcms_plate_(r)', 'exp_workflow_barcode', 'substrate_concentration_(r)'], plot_scatterplot=False, plot_histogram=False, calc_FIOP_bygrp=False, get_ctrls_from_vals=False):

        print('****************************')
        print('Computing Plate Analytics...')
        print('****************************')
        
        grp_list = sorted([l for l in list(set(df[groupby[0]])) if isinstance(l, str)])
        stats_table_dict, CTRL_stats, metricset_plots = {}, {}, {}
        CTRL_samples = {'pos':{}, self.neg_ctrltype:{}}
        for metric in self.metric_list:
            df.loc[:, metric] = df[metric].fillna(value=np.nan)

        ########################################################################
        # Get Analytic Metrics by MetricSet (e.g. Maldi, LcmsC18, Lcms Chiral) #
        ########################################################################
        
        for metricsetname, metricset in self.metric_dict.items():
            print(metricsetname, metricset)
            metricset_analysis = self.analytics_classes[metricsetname]
            stats_table_bymetricset, stats_bymetricset_CTRL_GRP, stats_bymetricset_CTRL_ALL, CTRL_samples_bymetricset = metricset_analysis.get_group_analytics(df, groupby, grp_list, metricset, calc_FIOP_bygrp=calc_FIOP_bygrp, get_ctrls_from_vals=get_ctrls_from_vals)
            
            # update joint containers
            stats_table_dict.update({metricsetname:stats_table_bymetricset})
            CTRL_stats.update(stats_bymetricset_CTRL_GRP)
            CTRL_samples['pos'].update(CTRL_samples_bymetricset['pos'])
            CTRL_samples[self.neg_ctrltype].update(CTRL_samples_bymetricset[self.neg_ctrltype])
            
            ## get plots by metricset ##
            metricset_plots[metricsetname] = metricset_analysis.get_plate_visualizations(df, metricset, groupby[0])

        ################################
        ## Compile Metric Stats Table ##
        ################################
        cols_to_get = groupby + ['n_var'] + [f'n_{metricset}' for metricset in self.metric_dict] + [f'{metric}{suffix}' for metric in self.metric_list for suffix in self.suffix_list]
        cols_to_display = groupby + ['n_var'] + [f'n_{metricset}' for metricset in self.metric_dict] + [f'{metric}{suffix}' for metric in self.metric_list for suffix in self.suffix_list_to_display]
        stats_table = self.get_analytics_table(stats_table_dict, cols_to_get, groupbyname='plate', sort_by=self.sort_by, cols_to_display=cols_to_display)
        
        #################################
        ## get plots across metricsets ##
        #################################
        metricset, metricnames = self.scatterplot_axes['EE_vs_Conversion']
        scatterplot, histogram = None, None
        # scatter plot
        if plot_scatterplot:
            scatterplot=plot_scatter_combined_metrics(df, metricset, metricnames, groupby[0], grp_list, groupbyname='plate', save_plot_to_s3=True, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder, fname_suffix=self.fname_suffix)
        # histogram
        if plot_histogram:
            histogram = plot_histogram_combined_metrics(df, CTRL_samples, metricset, metricnames, groupby[0], grp_list, groupbyname='plate', save_plot_to_s3=True, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder, fname_suffix=self.fname_suffix)
        
        # #########################
        # ## Create Panel Column ##
        # #########################
        analytics_panel = []
        if self.get_dashboard_panel:
            text = 'PLATE ANALYTICS'
            analytics_panel = self.create_panel_columns(text=text, metricset_plots=metricset_plots, scatterplot=scatterplot, histogram=histogram, stats_table=stats_table, groupby='plate', cols_to_display=cols_to_display)
        
        print('Obtained plate analytics metrics. \n')
        return stats_table, analytics_panel
    
    
    def GET_LIBRARY_ANALYTICS(self, df, display_table=True, groupby=['library_ref', 'library_barcode', 'exp_workflow_barcode'], plot_scatterplot=False, plot_histogram=True, calc_FIOP_bygrp=False, get_ctrls_from_vals=False):

        print('******************************')
        print('Computing Library Analytics...')
        print('******************************')
        
        grp_list = sorted([l for l in list(set(df[groupby[0]])) if isinstance(l, str)])
        grp_list_negctrl = [l for l in grp_list if any(substr in l.lower() for substr in ['negative', 'ev', 'empty'])]
        grp_list_posctrl = [l for l in grp_list if any(substr in l.lower() for substr in ['positive', 'wt'])]
        grp_list_expctrl = [l for l in grp_list if (l not in grp_list_negctrl) and (l not in grp_list_posctrl)]
        stats_table_dict, CTRL_stats, metricset_plots = {}, {}, {}
        CTRL_samples = {'pos':{}, self.neg_ctrltype:{}}
        for metric in self.metric_list:
            df.loc[:, metric] = df[metric].fillna(value=np.nan)

        ########################################################################
        # Get Analytic Metrics by MetricSet (e.g. Maldi, LcmsC18, Lcms Chiral) #
        ########################################################################
        
        for metricsetname, metricset in self.metric_dict.items():
            print(metricsetname, metricset)
            metricset_analysis = self.analytics_classes[metricsetname]
            stats_table_bymetricset, stats_bymetricset_CTRL_GRP, stats_bymetricset_CTRL_ALL, CTRL_samples_bymetricset = metricset_analysis.get_group_analytics(df, groupby,
            grp_list, metricset, calc_FIOP_bygrp=calc_FIOP_bygrp, get_ctrls_from_vals=get_ctrls_from_vals,
            filter_posctrl=('library_ref', grp_list_posctrl), filter_negctrl=('library_ref', grp_list_negctrl), filter_expctrl=('library_ref', grp_list_expctrl))
            
            # update joint containers
            stats_table_dict.update({metricsetname:stats_table_bymetricset})
            CTRL_stats.update(stats_bymetricset_CTRL_GRP)
            CTRL_samples['pos'].update(CTRL_samples_bymetricset['pos'])
            CTRL_samples[self.neg_ctrltype].update(CTRL_samples_bymetricset[self.neg_ctrltype])
            
            ## get plots by metricset ##
            metricset_plots[metricsetname] = metricset_analysis.get_library_visualizations(df, stats_table_bymetricset, metricset, groupby[1], 100)

        ################################
        ## Compile Metric Stats Table ##
        ################################
        cols_to_get = groupby + ['n_var'] + [f'n_{metricset}' for metricset in self.metric_dict] + [f'{metric}{suffix}' for metric in self.metric_list for suffix in self.suffix_list]
        cols_to_display = groupby + ['n_var'] + [f'n_{metricset}' for metricset in self.metric_dict] + [f'{metric}{suffix}' for metric in self.metric_list for suffix in self.suffix_list_to_display]
        stats_table = self.get_analytics_table(stats_table_dict, cols_to_get, groupbyname='library', sort_by=self.sort_by, cols_to_display=cols_to_display)
        
        #################################
        ## get plots across metricsets ##
        #################################
        metricset, metricnames = self.scatterplot_axes['EE_vs_Conversion']
        scatterplot, histogram = None, None
        # scatter plot
        if plot_scatterplot:
            scatterplot = plot_scatter_combined_metrics(df, metricset, metricnames, groupby[0], grp_list, groupbyname='library', save_plot_to_s3=True, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder, fname_suffix=self.fname_suffix)
        # histogram
        if plot_histogram:
            histogram = plot_histogram_combined_metrics(stats_table.iloc[3:], CTRL_samples, metricset, metricnames, groupby=groupby[0], groupbyname='library', plot_exp_metric_stats=['median', 'iqr'], showplot=False, neg_ctrltype='EV', fig=None, ax=None, save_plot_to_s3=True, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder, fname_suffix=self.fname_suffix)
        
        # #########################
        # ## Create Panel Column ##
        # #########################
        analytics_panel = []
        if self.get_dashboard_panel:
            text = 'LIBRARY ANALYTICS'
            analytics_panel = self.create_panel_columns(text=text, metricset_plots=metricset_plots, scatterplot=scatterplot, histogram=histogram, stats_table=stats_table, groupby='library', cols_to_display=cols_to_display)
        
        print('Obtained library analytics metrics. \n')
        return stats_table, analytics_panel
    
    
    def GET_UNIT_ANALYTICS(self, df, display_table=True, groupby=['exp_workflow_barcode', 'substrate_concentration_(r)', 'lcms_plate_(r)'], plot_scatterplot=True, plot_histogram=True, calc_FIOP_bygrp=False, get_ctrls_from_vals=True):
        
        print('***************************')
        print('Computing Unit Analytics...')
        print('***************************')
        
        # group list
        grp_list = sorted(list(set(df[groupby[0]])))
        stats_table_dict, CTRL_stats, metricset_plots = {}, {}, {}
        CTRL_samples = {'pos':{}, self.neg_ctrltype:{}}

        ########################################################################
        # Get Analytic Metrics by MetricSet (e.g. Maldi, LcmsC18, Lcms Chiral) #
        ########################################################################
        
        for metricsetname, metricset in self.metric_dict.items():
            print(metricsetname, metricset)
            metricset_analysis = self.analytics_classes[metricsetname]
            stats_table_bymetricset, stats_bymetricset_CTRL_GRP, stats_bymetricset_CTRL_ALL, CTRL_samples_bymetricset = metricset_analysis.get_group_analytics(df, groupby, grp_list, metricset, filter_posctrl=('ctrl_type', ['pos']), filter_negctrl=('ctrl_type', ['EV']), filter_expctrl=('ctrl_type', ['exp']), calc_FIOP_bygrp=calc_FIOP_bygrp, get_ctrls_from_vals=get_ctrls_from_vals)
            
            # update joint containers
            stats_table_dict.update({metricsetname:stats_table_bymetricset})
            CTRL_stats.update(stats_bymetricset_CTRL_GRP)
            CTRL_samples['pos'].update(CTRL_samples_bymetricset['pos'])
            CTRL_samples[self.neg_ctrltype].update(CTRL_samples_bymetricset[self.neg_ctrltype])
            
            ## get plots by metricset ##
            metricset_plots[metricsetname] = metricset_analysis.get_unit_visualizations(df, metricset, groupby[0])

        ################################
        ## Compile Metric Stats Table ##
        ################################
        cols_to_get = groupby + ['n_var'] + [f'n_{metricset}' for metricset in self.metric_dict] + [f'{metric}{suffix}' for metric in self.metric_list for suffix in self.suffix_list]
        cols_to_display = groupby + ['n_var'] + [f'n_{metricset}' for metricset in self.metric_dict] + [f'{metric}{suffix}' for metric in self.metric_list for suffix in self.suffix_list_to_display]
        stats_table = self.get_analytics_table(stats_table_dict, cols_to_get, groupbyname='unit', sort_by=self.sort_by, cols_to_display=cols_to_display)
        
        #################################
        ## get plots across metricsets ##
        #################################
        metricset, metricnames = self.scatterplot_axes['EE_vs_Conversion']
        scatterplot, histogram = None, None
        # scatter plot
        if plot_scatterplot:
            scatterplot = plot_scatter_combined_metrics(df, metricset, metricnames, groupby[0], grp_list, groupbyname='unit', save_plot_to_s3=True, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder, fname_suffix=self.fname_suffix)
        # histogram
        if plot_histogram:
            histogram = plot_histogram_combined_metrics(df, CTRL_samples, metricset, metricnames, groupby[0], grp_list, groupbyname='unit', save_plot_to_s3=True, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder, fname_suffix=self.fname_suffix)
        
        # #########################
        # ## Create Panel Column ##
        # #########################
        analytics_panel = []
        if self.get_dashboard_panel:
            text = 'UNIT ANALYTICS'
            analytics_panel = self.create_panel_columns(text=text, metricset_plots=metricset_plots, scatterplot=scatterplot, histogram=histogram, stats_table=stats_table, groupby='unit', cols_to_display=cols_to_display)
        
        print('Obtained unit analytics metrics. \n')
        return stats_table, analytics_panel

    
    def GET_EXPERIMENT_ANALYTICS(self, df, exp_column='exp_condition', groupby=['exp_condition', 'enzyme_barcode', 'mutations', 'library_ref', 'substrate_concentration_(r)', 'lcms_plate_(r)'], sort_by=None, display_table=True, plot_scatterplot=False, plot_histogram=False, calc_FIOP_bygrp=False, get_ctrls_from_vals=False, table_suffix=''):
        
        print('*********************************')
        print('Computing Experiment Analytics...')
        print('*********************************')
        
        if exp_column not in groupby: 
            groupby = [exp_column] + groupby
        grp_list = sorted([l for l in list(set(df[exp_column])) if isinstance(l, str)])
        stats_table_dict, CTRL_stats = {}, {}
        CTRL_samples = {'pos':{}, self.neg_ctrltype:{}}
        for metric in self.metric_list:
            df.loc[:, metric] = df[metric].fillna(value=np.nan)

        ########################################################################
        # Get Analytic Metrics by MetricSet (e.g. Maldi, LcmsC18, Lcms Chiral) #
        ########################################################################
        
        for metricsetname, metricset in self.metric_dict.items():
            print(metricsetname, metricset)
            metricset_analysis = self.analytics_classes[metricsetname]
            stats_table_bymetricset, stats_bymetricset_CTRL_GRP, stats_bymetricset_CTRL_ALL, CTRL_samples_bymetricset = metricset_analysis.get_group_analytics(df, groupby, grp_list, metricset, calc_FIOP_bygrp=calc_FIOP_bygrp, get_ctrls_from_vals=get_ctrls_from_vals)
            
            # update joint containers
            stats_table_dict.update({metricsetname:stats_table_bymetricset})
            CTRL_stats.update(stats_bymetricset_CTRL_GRP)
            CTRL_samples['pos'].update(CTRL_samples_bymetricset['pos'])
            CTRL_samples[self.neg_ctrltype].update(CTRL_samples_bymetricset[self.neg_ctrltype])

        ################################
        ## Compile Metric Stats Table ##
        ################################
        cols_to_get = groupby + [f'n_{metricset}' for metricset in self.metric_dict] + [f'{metric}{suffix}' for metric in self.metric_list for suffix in self.suffix_list]
        cols_to_display = groupby + [f'n_{metricset}' for metricset in self.metric_dict] + [f'{metric}{suffix}' for metric in self.metric_list for suffix in self.suffix_list_to_display]
        stats_table = self.get_analytics_table(stats_table_dict, cols_to_get, groupbyname=exp_column, sort_by=sort_by, table_suffix=table_suffix, cols_to_display=cols_to_display)
            
        #################################
        ## get plots across metricsets ##
        #################################
        metricset, metricnames = self.scatterplot_axes['EE_vs_Conversion']
        scatterplot, histogram = None, None
        # scatter plot
        if plot_scatterplot:
            scatterplot = plot_scatter_combined_metrics(df, metricset, metricnames, groupby[0], grp_list, groupbyname='library', save_plot_to_s3=True, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder, fname_suffix=self.fname_suffix)
        # histogram
        if plot_histogram:
            histogram = plot_histogram_combined_metrics(stats_table.iloc[3:], CTRL_samples, metricset, metricnames, groupby=groupby[0], groupbyname='enzyme', plot_exp_metric_stats=['median', 'iqr'], showplot=False, neg_ctrltype='EV', fig=None, ax=None, save_plot_to_s3=True, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder, fname_suffix=self.fname_suffix)
        
        # #########################
        # ## Create Panel Column ##
        # #########################
        analytics_panel = []
        if self.get_dashboard_panel:
            text = 'EXPERIMENT ANALYTICS'
            analytics_panel = self.create_panel_columns(text=text, scatterplot=scatterplot, histogram=histogram, stats_table=stats_table, groupby='variant', cols_to_display=cols_to_display)
            
        print('Obtained variant experiment metrics. \n')
        return stats_table, analytics_panel
    
    
    
    def get_activity_correlation_boxplots(self, df_reps, metric_list=['measured_nonbinary_score_(r)', 'pellet_OD', 'measured_nonbinary_sum_(r)'], metricname_list=['RacemicProduct', 'PelletOD', 'RacemicSum'], groupby=['mutations'], var_col_idx=0, var_list=None, fname_prefix='', table_suffix=''):
        
        # get var_list to ORDER boxplot from dataframe order
        if var_list is None:
            var_list = [] 
            for var in df_reps[groupby[var_col_idx]].tolist():
                if var not in var_list: 
                    var_list.append(var)

        # get groups
        df_orderedgrps_dict = {metric: OrderedDict() for metric in metric_list}
        
        # loop over individual variants in var_list
        for var in var_list: 
            df_var = df_reps.loc[df_reps[groupby[var_col_idx]]==var, groupby+metric_list]
            
            # APPLY ADDITIONAL FILTERS TO GET SUBGROUPS (e.g. EXP_CONDITION) 
            if len(groupby)>1:
                grp_list = df_var[groupby].drop_duplicates().sort_values(by=groupby).to_records(index=False).tolist()
                for grp in grp_list: 
                    df_grp = df_var.copy()
                    for col, val in zip(groupby, grp): 
                        df_grp = df_grp.loc[df_grp[col]==val]
                    
                    # update dict for that metric
                    for metric in metric_list:
                        df_orderedgrps_dict[metric][grp] = df_grp[metric]
                        
            # NO ADDITIONAL FILTERS
            ## update dict for that metric
            else:
                for metric in metric_list:
                    df_orderedgrps_dict[metric][var] = df_var[metric]
    
        # initialize plot list
        plot_list = []
        
        # plot boxplot for each metric
        for i, (metric, metricname) in enumerate(zip(metric_list, metricname_list)): 
            
            df_orderedgrps = pd.DataFrame(df_orderedgrps_dict[metric])
            boxplot_list = plot_boxplot(df_orderedgrps, [metric], [metricname], groupby=None, fname_prefix=fname_prefix, fname_suffix=f'_{"-".join(groupby)}', img_format='png', plate_barcode=None, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder, save_plot_to_s3=False, save_plot_to_lims=False, show_n=True)
            plot_list += boxplot_list
            
            # save plots to S3
            if self.save_plots_to_s3:
                img_fname = f'BOXPLOT_{fname_prefix}{metricname}_{"-".join(groupby)}'
                img_header = f'{fname_prefix}{metricname} boxplots grouped by {", ".join(groupby)}'
                self.save_content_to_s3(img_fname, file_format='png', content_header=img_header)
        
        return plot_list
    
    
    def get_activity_correlation_scatterplots(self, df_reps, xmetric_list=['pellet_OD', 'measured_nonbinary_sum_(r)'], xmetricname_list=['PelletOD', 'RacemicSum'], ymetric='measured_nonbinary_score_(r)', ymetricname='RacemicProduct', groupby='mutations', table_suffix='', showlegend=False):
        
        # filter stats table to remove null entries for column to groupby
        df_reps = df_reps.loc[~df_reps[groupby].isnull()]
        
        # initialize plot list
        plot_list = []
        
        # get scatterplots
        for xmetric, xmetricname in zip(xmetric_list, xmetricname_list):
            if xmetric in df_reps:
                img_fname = f'SCATTERPLOT_{ymetricname}-vs-{xmetricname}_{groupby}'
                scatterplot = plot_scatter_combined_metrics(df_reps, metricset=[xmetric, ymetric], metricnames=[xmetricname, ymetricname], groupby=groupby, save_plot_to_s3=False, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder, fname_suffix=f'{self.fname_suffix}{table_suffix}', aspect_ratio=1, showlegend=showlegend)
                plot_list.append(scatterplot)
                
                # save plots to S3
                if self.save_plots_to_s3:
                    img_fname = f'SCATTERPLOT_{ymetricname}-vs-{xmetricname}_{groupby}{self.fname_suffix}{table_suffix}'
                    img_header = f'{ymetricname} vs. {xmetricname} scatterplot grouped by {groupby}'
                    self.save_content_to_s3(img_fname, file_format='png', content_header=img_header)
                
        return plot_list
    
        
    def get_activity_CV_correlation_scatterplots(self, stats_table_filt, xmetric_list=['pellet_OD_cv', 'measured_nonbinary_sum_(r)_cv'], xmetricname_list=['PelletOD-CV', 'RacemicSum-CV'], ymetric='measured_nonbinary_score_(r)_cv', ymetricname='RacemicProduct-CV', groupby='mutations',table_suffix=''):

        # filter stats table
        stats_table_filt = stats_table_filt.loc[~stats_table_filt[groupby].isnull()]
        
        # initialize plot list
        plot_list = []
        
        for xmetric, xmetricname in zip(xmetric_list, xmetricname_list):
            if xmetric in stats_table_filt:
                # scatterplot
                scatterplot = plot_scatter_combined_metrics(stats_table_filt, metricset=[xmetric, ymetric], metricnames=[xmetricname, ymetricname], groupby=groupby, save_plot_to_s3=True, s3_bucket=self.s3_bucket, s3_subfolder=self.s3_subfolder, fname_suffix=f'{self.fname_suffix}{table_suffix}', aspect_ratio=1, showlegend=False)
                
            else:
                scatterplot = None
                print(f'{xmetric} not found in dataframe.')
            plot_list.append(scatterplot)
        return plot_list

    
    def get_analytics_table(self, stats_table_dict, cols_to_get, groupbyname='unit',
                            display_table=True, include_ctrl_stats=True, sort_by=None, table_suffix='', cols_to_display=None):
        """
        Get statistics of analytics metrics in table form
        """
        metricsetname_list = list(self.metric_dict.keys())
        for i, metricsetname in enumerate(metricsetname_list):
            if i==0:
                # stats_table = pd.DataFrame(stats_table_dict[metricsetname])
                stats_table = stats_table_dict[metricsetname]
            else:
                # stats_table_toconcat = pd.DataFrame(stats_table_dict[metricsetname])
                stats_table_toconcat = stats_table_dict[metricsetname]
                stats_table = pd.concat([stats_table, stats_table_toconcat], axis=1)
                
        # remove duplicate columns & drop empty columns
        if len(metricsetname_list) > 1:
            stats_table = stats_table[[c for c in cols_to_get if c in stats_table]]
        stats_table = stats_table.loc[:,~stats_table.columns.duplicated()].dropna(axis=1, how='all')
        # sort exp rows of table
        if sort_by is not None:
            metric, stat = sort_by[0], sort_by[1]
            stats_table = pd.concat([stats_table.iloc[:3], stats_table.iloc[3:].sort_values(f'{metric}_{stat}', ascending=False).reset_index(drop=True)])
        
        # columns to display that are in table
        cols_to_display = [c for c in cols_to_display if c in stats_table]
        
        if display_table and not self.get_dashboard_panel:
            pd.set_option('display.float_format', '{:.3f}'.format)
            display(stats_table[cols_to_display])
            
        # save table to S3
        csv_fname = f'{groupbyname}_METRICS{self.fname_suffix}{table_suffix}'
        csv_header = f'Table of metric statistics (grouped by {groupbyname})'
        self.save_content_to_s3(csv_fname, stats_table, file_format='csv', content_header=csv_header, content_type='dataframe', cols_to_display=cols_to_display)

        return stats_table
    
    
    def update_analytics_table_postgres(self, df, enzyme_analytics_table='combi_analytics_table', update_id=True, update_analytics_table=False):
        
        # update id column
        if update_id:
            ids = [f'{w}_{p}_{a}' for w, p, a in zip(df['exp_workflow_barcode'].tolist(), df['lcms_plate_(r)'].tolist(), df['lcms_address_(r)'].tolist())]
            df['id'] = ids
        
        # add missing columns to table
        df = DatabaseInterface(table=enzyme_analytics_table).pad_with_missing_columns(df, combi_analytics_metrics)
        # remove rows which don't have any associated source plates
        df = df[~df['source_plate_(+)'].isnull() | ~df['source_plate_(-)'].isnull() | ~df['source_plate_(r)'].isnull()]
        
        # update enzyme_analytics_table records in S3 by unique source address, source plate combinations
        if update_analytics_table:
            source_plate_list = list(set(df['source_plate_(r)'].tolist()))
            print('source_plate_(r) list:', source_plate_list)
            for source_plate in source_plate_list:
                # get data relevant to source plate
                df_plate = df[df['source_plate_(r)']==source_plate].copy()
                # clear previous data for source plate
                db_detections = DatabaseInterface(table=enzyme_analytics_table)
                db_detections.clear_matching_records({'source_plate_(r)': source_plate})
                # update with new data
                db_detections.table_add(df_plate)
                print(f'Updated postgres enzyme_analytics_table for all (+), (-), (r) samples associated with the variants on source plate {source_plate}.')
        return df
    
    
    def save_content_to_s3(self, fname, content=None, file_format='png', content_header=None, content_type='image', cols_to_display=None):
        filename = f'{fname}.{file_format}'
        if content_type=='image':
            s3_imgupload(self.s3_bucket, f'{self.s3_subfolder}{filename}', dpi=300)
            self.dashboard_contents.append({'filename':filename, 'text':content_header, 'content_type':content_type})
        elif content_type=='dataframe':
            s3_df2csv(content, self.s3_bucket, f'{self.s3_subfolder}{filename}')
            self.dashboard_contents.append({'filename':filename, 'text':content_header, 'content_type':content_type, 'cols_to_display':cols_to_display})
            
        print(f"Saved {self.s3_subfolder}{filename} to S3.")
    
    
    def fetch_saved_panel_data(self, s3_bucket=None, s3_subfolder=None, dashboard_contents=None, fig_mindim=5, dpi=150, cols_to_display=None):
        """
        Combine plots and tables into a Dashboard Panel
        """
        from io import StringIO
            
        # fetch dashboard contents
        if s3_bucket is None:
            s3_bucket = self.s3_bucket
        if s3_subfolder is None:
            s3_subfolder = self.s3_subfolder
        if dashboard_contents is None: 
            filename = f'{s3_subfolder}dashboard_contents.pkl'
            dashboard_contents = download_from_s3(s3_bucket, filename)
            
        print('Re-loading pre-calculated data...')
        panel_columns = []

        for item in dashboard_contents:
            
            # filename
            filename = item['filename']
            fname = f"{s3_subfolder}{filename}"
            
            # text
            text = item['text'] if 'text' in item else None
            if text is not None:
                panel_columns.append(pn.Column(f'#{text}'))
            
            # image
            if item['content_type']=='image':
                fig = None
                try:
                    figdata = s3_imgdownload(s3_bucket, fname)
                    figdata_dim = figdata.shape
                    figdata_mindim = min(figdata_dim[0], figdata_dim[1])
                    fig_mindim = figdata_mindim/300 
                    aspect_ratio = figdata_dim[0] / figdata_dim[1]
                    if aspect_ratio < 1: # wide
                        figsize = (fig_mindim/aspect_ratio, fig_mindim)
                    else:
                        figsize = (fig_mindim, aspect_ratio*fig_mindim)
                    fig = plt.figure(figsize=figsize, dpi=dpi)
                    plt.imshow(figdata)
                    plt.axis('off')
                except Exception as e:
                    print(e, f'\n Could not load {filename} from S3.')
                else:
                    # create panel
                    panel_columns.append(pn.Column(pn.pane.Matplotlib(fig), scroll=True, background="White"))
                    print(f'Loaded {filename} from S3.')     
            
            # dataframe
            elif item['content_type']=='dataframe':
                try:
                    stats_table = s3_csv2df(s3_bucket, fname, skiprows=None, index_col=0, header=0)
                    df_widget_statstable = pn.widgets.DataFrame(stats_table.round(4), name='DataFrame', autosize_mode='fit_viewport')
                    cols_to_display = item['cols_to_display']
                    # create widget for download
                    sio_statstable = StringIO()
                    stats_table.to_csv(sio_statstable)
                    sio_statstable.seek(0)
                    csv_download_statstable = pn.widgets.FileDownload(file=sio_statstable, filename=filename, button_type="primary")
                    # table
                    if cols_to_display is not None:
                        stats_table = stats_table[[c for c in cols_to_display if c in stats_table]]
                except Exception as e:
                    print(e, f'\n Could not load {filename} from S3.')     
                else:
                    panel_columns.append(pn.Column(csv_download_statstable)) 
                    panel_columns.append(pn.Column(df_widget_statstable, height=600, scroll=True, background="White"))
                    print(f'Loaded {filename} from S3.')     
        
        return panel_columns
    
    
    def create_panel_columns(self, text=None, metricset_plots={}, scatterplot=None, histogram=None, stats_table=None, data_reps=None, groupby='ctrl_type', cols_to_display=None):
        """
        Combine plots and tables into a Dashboard Panel
        """
        from io import StringIO
        
        panel_columns = []
        
        # text
        if text is not None:
            heading = pn.pane.Markdown(f'# {text}', style={'color': 'White', 'background-color':'Black', 'font-size':'32-pt'})
            panel_columns.append(pn.Column(heading, background='Black'))
        
        # metricset plots
        for metricset, metricset_plotlist in metricset_plots.items():
            metricset_plot_pane = [pn.pane.Matplotlib(fig) for fig in metricset_plotlist]
            panel_columns.append(pn.Column(f'# {metricset} metric plots', *metricset_plot_pane, scroll=True, background="White"))

        # scatter plot + histogram
        if scatterplot is not None or histogram is not None:
            scatterplot_pane = pn.pane.Matplotlib(scatterplot)
            histogram_pane = pn.pane.Matplotlib(histogram)
            panel_columns.append(pn.Column('# Scatterplot & Histogram of Conversion & Enantioselectivity', pn.Row(scatterplot_pane, histogram_pane, scroll=True, background="White")))
        
        # stats table
        if stats_table is not None:
            sio_statstable = StringIO()
            stats_table.to_csv(sio_statstable)
            sio_statstable.seek(0)
            # download button
            csv_download_statstable = pn.widgets.FileDownload(file=sio_statstable, filename=f'{groupby}_metricstats.csv', button_type="primary")
            # table
            if cols_to_display is not None:
                stats_table = stats_table[[c for c in cols_to_display if c in stats_table]]
            df_widget_statstable = pn.widgets.DataFrame(stats_table.round(4), name='DataFrame', autosize_mode='fit_viewport')
            panel_columns.append(pn.Column(f'# Table of metric statistics (grouped by {groupby})', csv_download_statstable)) 
            panel_columns.append(pn.Column(df_widget_statstable, height=600, scroll=True, background="White"))
        
        # data reps table
        if data_reps is not None:
            sio_datareps = StringIO()
            data_reps.to_csv(sio_datareps)
            sio_datareps.seek(0)
            # download button
            csv_download_datareps = pn.widgets.FileDownload(file=sio_datareps, filename=f'{groupby}_datareps.csv', button_type="primary")
            panel_columns.append(pn.Column(f'# Replicate data (grouped by {groupby})', csv_download_datareps, width=1600, scroll=False, background="White"))

        return panel_columns 