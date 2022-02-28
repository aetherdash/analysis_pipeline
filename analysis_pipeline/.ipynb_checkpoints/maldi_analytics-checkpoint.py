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

from analytics_utils.analytics_utils.table_properties import *
from analytics_utils.analytics_utils.get_analytics_metrics import get_ctrl_vals, get_lcms_qc_derivedmetrics_checks
from analytics_utils.analytics_utils.db_interface import DatabaseInterface
from analytics_utils.analytics_utils.s3_interface import download_from_s3, upload_to_s3, s3_imgupload, s3_df2csv
from analytics_utils.analytics_utils.visualization_utils import generate_plots, df2array_dict
from analytics_utils.analytics_utils.lims_utils import lims_post_matplotlib_img
from messages.detectionspipeline_producers import ML_producer

        
class MaldiAnalytics: 
    
    def __init__(self, 
                 model_bucket='ml-models-registry',
                 model_subfolder='exerevnetes-preprocessing-models/ProjectX_UnitX/',
                 model_fname='model',
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

        self.model_bucket = model_bucket
        self.model_subfolder = model_subfolder
        self.model_fname = model_fname
        self.fname_suffix = fname_suffix
        self.model_path = f'{self.model_subfolder}{self.model_fname}'
        self.local_dir = local_dir
        self.cols_to_get = {'maldi_detections_1': maldi_detections_columns_to_get,
                            'maldi_detections_2': maldi_detections_columns_to_get,
                            'lcms_detections': lcms_detections_columns_to_get,
                            'lcms_detections_chiral': lcms_detections_chiral_columns_to_get}
        self.combine_maldi_lcms_by = combine_maldi_lcms_by
        self.posctrl_enzyme_barcode = posctrl_enzyme_barcode
        self.plot_enantiomer_lcms_metrics = plot_enantiomer_lcms_metrics
        self.col_suffix_db = col_suffix_db 
        self.save_plots_to_lims = save_plots_to_lims
        self.get_2sample_from_1sample_dataset = get_2sample_from_1sample_dataset
        self.ctrl_to_color_mapping = {'exp':'r', 'pos':'g', 'EV':'k'}
    
    
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
    
    
#     def GET_PLATE_ANALYTICS:
#         for plate in self.plate_list:

#             # get maldi plate analytics
#             df_plate = data[self.maldi_detections_1][data[self.maldi_detections_1][f'plate{col_suffix}']==plate].copy()
#             if len(df_plate) > 0:

#                 ###############################
#                 #### maldi plate analytics ####
#                 ###############################
#                 if update_plate_analytics_table:
#                     maldi_plate_analytics_dict = get_maldi_plate_analytics(df_plate, neg_ctrltype=self.neg_ctrltype, col_suffix=substrate_barcode, 
#                                                                          use_predicted_labels=predict_labels)
#                     plate_analytics_df.append(maldi_plate_analytics_dict)
#                     plate_indices.append(plate)
#                     print(f'Obtained plate analytics dataframe for {plate}.')         
#                     # clear previous data for source plate
#                     db_detections = DatabaseInterface(table=plate_analytics_table)
#                     db_detections.clear_matching_records({'plate': plate})
#                     # update with new data
#                     db_detections.table_add(pd.DataFrame(maldi_plate_analytics_dict, index=[0])[maldi_plate_analytics_metrics])
#                     print(f'Updated postgres plate_analytics_table for {plate}.')     


#                 ##############################
#                 # plot boxplots for controls #
#                 ##############################
#                 if generate_boxplots:
#                     columns_to_plot=[f'ptr{col_suffix}', f'pk_prod_stdz_172{col_suffix}', f'pk_sub_stdz_172{col_suffix}']
#                     boxplot = self.plot_boxplot(df_plate, plate, columns_to_plot, 
#                                                 file_suffix='_maldi-controls', s3_bucket_img='maldi-file-store')

#                 ##########################
#                 ## plot plate heat maps ##
#                 ##########################
#                 if generate_heatmaps_measured:
#                     for col_suffix in self.split_by_list:

#                         substrate_barcode = col_suffix
#                         col_suffix = col_suffix if col_suffix=='' else f'_{col_suffix}'
#                         label_condition_suffix = '_maldi' if get_lcms_labels else ''
#                         col_binary = f'true_score_binary{col_suffix}{label_condition_suffix}'
#                         col_nonbinary = f'{self.lcms_output_feature}{col_suffix}'

#                         # generate visualizations and post to maldi-file-store    
#                         for plate in self.plate_list:
#                             if substrate_barcode=='' or substrate_barcode in self.substrate_list[plate]:
#                                 # BINARY (1-sample)
#                                 if self.maldi_detections_1 in data and col_binary in data[self.maldi_detections_1]:
#                                     df_plate = data[self.maldi_detections_1][data[self.maldi_detections_1][f'plate{col_suffix}']==plate].copy()
#                                     self.generate_visualizations(df_plate, 'maldi1sample_classification', 
#                                                                  cols_to_plot=[col_binary], 
#                                                                  spectra_ids_col=f'spectra_ids{col_suffix}',
#                                                                  plate_barcode=plate, file_suffix='-1_binary_measured')


#                                 # BINARY (2-sample)
#                                 if self.maldi_detections_2 in data and col_binary in data[self.maldi_detections_2]:
#                                     df_plate = data[self.maldi_detections_2][data[self.maldi_detections_2][f'plate{col_suffix}']==plate].copy()
#                                     self.generate_visualizations(df_plate, 'maldi2sample_classification', 
#                                                                  cols_to_plot=[col_binary], 
#                                                                  spectra_ids_col=f'spectra_ids{col_suffix}',
#                                                                  plate_barcode=plate, file_suffix='-2_binary_measured')
#                                     if get_lcms_labels: 
#                                         self.generate_visualizations(df_plate, 'maldi2sample_classification', 
#                                                                      cols_to_plot=[f'true_score_binary{col_suffix}_lcms'], 
#                                                                      spectra_ids_col=f'spectra_ids{col_suffix}',
#                                                                      plate_barcode=plate, file_suffix='-2_binary_measured_lcmslabeled')                                

#                                 # NON-BINARY
#                                 if 'maldi_lcms' in data and col_nonbinary in data['maldi_lcms']:
#                                     df_plate = data['maldi_lcms'][data['maldi_lcms'][f'plate{col_suffix}']==plate].copy()
#                                     self.generate_visualizations(df_plate, 'maldi2lcms_regression', 
#                                                              cols_to_plot=[col_nonbinary], 
#                                                              spectra_ids_col=f'spectra_ids{col_suffix}',
#                                                              plate_barcode=plate, file_suffix='_nonbinary_measured')


#                 # generate heatmaps of predicted labels
#                 if generate_heatmaps_predicted and enzyme_analytics_df is not None :
#                     if substrate_barcode=='' or substrate_barcode in self.substrate_list[plate]:
#                         # BINARY
#                         df_plate = enzyme_analytics_df[enzyme_analytics_df[f'maldi_plate{col_suffix_db}']==plate].copy()
#                         self.generate_visualizations(df_plate, 'maldi2sample_classification', 
#                                                      cols_to_plot=[f'predicted_binary_score{col_suffix_db}'], 
#                                                      spectra_ids_col=f'maldi_address{col_suffix_db}',
#                                                      plate_barcode=plate, file_suffix='_binary_pred')

#                         # NON-BINARY
#                         self.generate_visualizations(df_plate, 'maldi2lcms_regression', 
#                                                  cols_to_plot=[f'predicted_nonbinary_score{col_suffix_db}'], 
#                                                  spectra_ids_col=f'maldi_address{col_suffix_db}',
#                                                  plate_barcode=plate, file_suffix='_nonbinary_pred') 

                            
    

    def GET_ENZYME_ANALYTICS(self, enzyme_analytics_df, split_by_list, predict_labels=False):

        ###################################
        # generate enzyme analytics plots #
        ###################################    
          
        # plot enantiomer MEASURED lcms metrics
        if 'measured_nonbinary_score_(+)' in enzyme_analytics_df and self.plot_enantiomer_lcms_metrics is not None:
            plot_lcms_metrics_for_ctrls(enzyme_analytics_df, 
                                        x_suffix='(-)',
                                        y_suffix='(+)', 
                                        ctrls_to_plot=['exp','pos',self.neg_ctrltype], 
                                        metric_list=['measured_nonbinary_score'])
            # upload plot to S3
            s3_imgupload(self.model_bucket, f"{self.model_subfolder}lcms_enantiomer_metrics_MEASURED{self.fname_suffix}.jpg")
            plt.show()
        
        # plot enantiomer PREDICTED lcms metrics
        if len(split_by_list)>1 and self.plot_enantiomer_lcms_metrics is not None and predict_labels:
            plot_lcms_metrics_for_ctrls(enzyme_analytics_df, 
                                        x_suffix='(-)',
                                        y_suffix='(+)', 
                                        ctrls_to_plot=['exp','pos',self.neg_ctrltype], 
                                        metric_list=['predicted_nonbinary_score'])
            # upload plot to S3
            s3_imgupload(self.model_bucket, f"{self.model_subfolder}lcms_enantiomer_metrics_PRED{self.fname_suffix}.jpg")
            plt.show()
    
    def get_derived_metrics(self):
        return None
    
    
    def get_qc_metrics(self):
        return None
    
    
    def get_plate_analytics(self):
        return None

    
    def GET_GROUP_ANALYTICS(self, df, groupby): 
        return None
    
    
    def GET_GROUP_ANALYTICS_PLOTS(self, df, groupby):
        return None

            
    def generate_visualizations(self, df, model_type, cols_to_plot=[], spectra_ids_col='spectra_ids', plate_barcode=None,
                               s3_bucket_img='maldi-file-store', file_suffix=''):
        
        if model_type=='maldi1sample_classification':
            df_plot = df
            set_logscale = False
        if model_type=='maldi2sample_classification' or model_type=='maldi2lcms_regression':
            columns = list(df.columns)
            df_plot = pd.DataFrame(columns=columns)
            for i in range(len(df)): 
                dfrow = df.iloc[i:i+1]
                addresses = str(dfrow.iloc[0][spectra_ids_col]).split(', ')
                if spectra_ids_col == 'spectra_ids':
                    addresses = [a[:a.find('_')] for a in addresses]
                df_to_append = pd.DataFrame(np.repeat(dfrow.values, len(addresses), axis=0), columns=columns)
                df_to_append['address'] = addresses
                df_plot = df_plot.append(df_to_append)
        generate_plots(df2array_dict(df_plot, cols_to_plot), 
                       plate_barcode=plate_barcode, 
                       plottype=f'{model_type}{file_suffix}', 
                       set_logscale=False, title_fontsize=14, figsize_factor=0.7)
        
        # upload image to S3
        img_fname = f'{plate_barcode}{file_suffix}.jpg'
        s3_imgupload(s3_bucket_img, f"{plate_barcode}/{img_fname}")
        
        # post plot to LIMS file store
        if self.save_plots_to_lims:
            sub_url = 'plates/' + plate_barcode + '/media/'
            lims_post_matplotlib_img(img_fname, sub_url)
            print(f"Posted {img_fname} plot to LIMS via requests module.")
            
        plt.show()
        return df_plot

    
    def plot_boxplot(self, df, plate_barcode, columns_to_plot, file_suffix='', s3_bucket_img='maldi-file-store'):
        
        fig = plt.figure();
        boxplot = df.boxplot(column=columns_to_plot, by=['ctrl_type'], 
                             layout=(1,len(columns_to_plot)), figsize=(7*len(columns_to_plot),7))
        
        # upload image to S3
        img_fname = f'{plate_barcode}{file_suffix}_boxplot.jpg'
        s3_imgupload(s3_bucket_img, f"{plate_barcode}/{img_fname}")
        
        # post plot to LIMS file store
        if self.save_plots_to_lims:
            sub_url = 'plates/' + plate_barcode + '/media/'
            lims_post_matplotlib_img(img_fname, sub_url)
            print(f"Posted {img_fname} plot to LIMS via requests module.")
            
        plt.show()
        return boxplot
    
    
def get_maldi_plate_analytics(df, neg_ctrltype='EV', col_suffix='', 
                              true_binary_score_colname='true_score_binary', predicted_binary_score_colname='predicted_binary_score', use_predicted_labels=False): 
    
    """
    Get analysis metrics for a plate
    """
    # create copy of dataframe and remove col_suffix from colnames if applicable
    col_suffix = col_suffix if col_suffix=='' else f'_{col_suffix}'
    cols_to_rename = [c for c in list(df.columns) if c.find(col_suffix)>-1]
    df = df.rename(columns={c: c.replace(col_suffix, '') for c in cols_to_rename})
    
    # column names 
    plate_meta_cols = ['plate', 'run', 'dev_or_prod', 'proj_name', 'proj_barcode', 'exp_workflow_name', 'exp_workflow_barcode']
    mtx_cols = ['pk_mtx_172_raw', 'pk_mtx_190_raw', 'pk_mtx_379_raw']
    rxn_cols = ['pk_sub_stdz_379', 'pk_prod_stdz_379', 'ptr', 'tstat_pk_sub_stdz_379', 'tstat_pk_prod_stdz_379', 'tstat_ptr']

    # segment data into different control types
    df_segmented = {}
    df_segmented['mtx'] = df[df.ctrl_type=='mtx']
    df_segmented['neg'] = df[df.ctrl_type==neg_ctrltype]
    df_segmented['pos'] = df[df.ctrl_type=='pos']
    df_segmented['exp'] = df[df.ctrl_type=='exp']

    # initialize plate_analytics_dict
    plate_analytics_dict = {c: None for c in maldi_plate_analytics_metrics}
    for metadata in plate_meta_cols:
        plate_analytics_dict[metadata] = df.iloc[0][metadata]
    for ctrltype in ['mtx', 'neg', 'pos', 'exp']:
        plate_analytics_dict[f'num_{ctrltype}_ctrl'] = int(len(df_segmented[ctrltype]))
    plate_analytics_dict['num_exp_variants'] = int(len(list(set(df_segmented['exp']['enzyme_barcode']))))

    # matrix controls
    for col in mtx_cols:
        df_subset = df_segmented['mtx']
        if len(df_subset) > 0:
            data = np.array(list(df_subset[col]))
            stats = calculate_stats_for_array(data)
            for metric in stats:
                plate_analytics_dict[f"{col}_{metric}"] = stats[metric]
    
    # negative, positive, experimental controls BINARY STATS (measured or predicted)
    if use_predicted_labels:
        binary_score_colname = predicted_binary_score_colname
    else: 
        binary_score_colname = true_binary_score_colname
    plate_analytics_dict['neg_ctrl_fpr'] = np.mean(np.array(list(df_segmented['neg'][binary_score_colname]))) if len(df_segmented['neg']) > 0 else np.nan
    plate_analytics_dict['pos_ctrl_fnr'] = 1-np.mean(np.array(list(df_segmented['pos'][binary_score_colname]))) if len(df_segmented['pos']) > 0 else np.nan
    plate_analytics_dict['exp_ctrl_hit_rate'] = np.mean(np.array(list(df_segmented['exp'][binary_score_colname]))) if len(df_segmented['exp']) > 0 else np.nan

    # negative, positive, experimental controls PEAK STATS
    for ctrltype in ['neg', 'pos', 'exp']:
        df_subset = df_segmented[ctrltype]
        if len(df_subset) > 0:
            plate_analytics_dict[f'{ctrltype}_ctrl_enzyme_barcode'] = ','.join(list(set(df_subset.enzyme_barcode.astype(str))))
            plate_analytics_dict[f'{ctrltype}_ctrl_substrate_barcode'] = ','.join(list(set(df_subset.substrate_barcode.astype(str))))
            for col in rxn_cols:
                data = np.array(list(df_subset[col]))
                stats = calculate_stats_for_array(data)
                for metric in stats:
                    plate_analytics_dict[f"{ctrltype}_ctrl_{col.replace('pk_','')}_{metric}"] = stats[metric]    

    return plate_analytics_dict