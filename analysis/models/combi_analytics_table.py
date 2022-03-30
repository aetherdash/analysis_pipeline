from django.db import models

class CombiAnalyticsTable(models.Model):

    class Meta:
        db_table = 'combi_analytics_table'
        
#    id = models.TextField()
    exp_workflow_barcode = models.TextField()
    exp_workflow_name = models.TextField()
    proj_barcode = models.TextField()
    proj_name = models.TextField()
    maldi_run = models.TextField()
    lcms_C18_run = models.TextField()
    lcms_chiral_run = models.TextField()
    maldi_plate_racemic = models.TextField(db_column='maldi_plate_(r)')
    maldi_address_racemic = models.TextField(db_column='maldi_address_(r)')
    lcms_plate_racemic = models.TextField(db_column='lcms_plate_(r)')
    lcms_address_racemic = models.TextField(db_column='lcms_address_(r)')
    source_plate_racemic = models.TextField(db_column='source_plate_(r)')
    source_address_racemic = models.TextField(db_column='source_address_(r)')
    maldi_plate_plus = models.TextField(db_column='maldi_plate_(+)')
    maldi_address_plus = models.TextField(db_column='maldi_address_(+)')
    lcms_plate_plus = models.TextField(db_column='lcms_plate_(+)')
    lcms_address_plus = models.TextField(db_column='lcms_address_(+)')
    source_plate_plus = models.TextField(db_column='source_plate_(+)')
    source_address_plus = models.TextField(db_column='source_address_(+)')
    maldi_plate_minus = models.TextField(db_column='maldi_plate_(-)')
    maldi_address_minus = models.TextField(db_column='maldi_address_(-)')
    lcms_plate_minus = models.TextField(db_column='lcms_plate_(-)')
    lcms_address_minus = models.TextField(db_column='lcms_address_(-)')
    source_plate_minus = models.TextField(db_column='source_plate_(-)')
    source_address_minus = models.TextField(db_column='source_address_(-)')
    ctrl_type = models.TextField()
    exp_condition = models.TextField()
    enzyme_barcode = models.TextField()
    sequence = models.TextField()
    mutations = models.TextField()
    hamming = models.IntegerField()
    reference_enzyme = models.TextField()
    substrate_barcode_racemic = models.TextField(db_column='substrate_barcode_(r)')
    substrate_barcode_plus = models.TextField(db_column='substrate_barcode_(+)')
    substrate_barcode_minus = models.TextField(db_column='substrate_barcode_(-)')
    substrate_smiles_racemic = models.TextField(db_column='substrate_smiles_(r)')
    substrate_smiles_plus = models.TextField(db_column='substrate_smiles_(+)')
    substrate_smiles_minus = models.TextField(db_column='substrate_smiles_(-)')
    substrate_concentration_racemic = models.FloatField(db_column='substrate_concentration_(r)')
    substrate_concentration_plus = models.FloatField(db_column='substrate_concentration_(+)')
    substrate_concentration_minus = models.FloatField(db_column='substrate_concentration_(-)')
    predicted_binary_score_racemic = models.FloatField(db_column='predicted_binary_score_(r)')
    predicted_nonbinary_score_racemic = models.FloatField(db_column='predicted_nonbinary_score_(r)')
    measured_nonbinary_score_racemic = models.FloatField(db_column='measured_nonbinary_score_(r)')
    measured_nonbinary_sum_racemic = models.FloatField(db_column='measured_nonbinary_sum_(r)')
    measured_conversion_racemic = models.FloatField(db_column='measured_conversion_(r)')
    predicted_binary_score_variance_racemic = models.FloatField(db_column='predicted_binary_score_variance_(r)')
    predicted_nonbinary_score_variance_racemic = models.FloatField(db_column='predicted_nonbinary_score_variance_(r)')
    predicted_binary_score_plus = models.FloatField(db_column='predicted_binary_score_(+)')
    predicted_nonbinary_score_plus = models.FloatField(db_column='predicted_nonbinary_score_(+)')
    measured_nonbinary_score_plus = models.FloatField(db_column='measured_nonbinary_score_(+)')
    measured_conversion_plus = models.FloatField(db_column='measured_conversion_(+)')
    predicted_binary_score_variance_plus = models.FloatField(db_column='predicted_binary_score_variance_(+)')
    predicted_nonbinary_score_variance_plus = models.FloatField(db_column='predicted_nonbinary_score_variance_(+)')
    predicted_binary_score_minus = models.FloatField(db_column='predicted_binary_score_(-)')
    predicted_nonbinary_score_minus = models.FloatField(db_column='predicted_nonbinary_score_(-)')
    measured_nonbinary_score_minus = models.FloatField(db_column='measured_nonbinary_score_(-)')
    measured_conversion_minus = models.FloatField(db_column='measured_conversion_(-)')
    predicted_binary_score_variance_minus = models.FloatField(db_column='predicted_binary_score_variance_(-)')
    predicted_nonbinary_score_variance_minus = models.FloatField(db_column='predicted_nonbinary_score_variance_(-)')
    measured_enantiomeric_excess_plus_over_minus = models.FloatField(db_column='measured_enantiomeric_excess_(+over-)')
    predicted_enantiomeric_excess_plus_over_minus = models.FloatField(db_column='predicted_enantiomeric_excess_(+over-)')
    measured_enantiomeric_ratio_plus_over_minus = models.FloatField(db_column='measured_enantiomeric_ratio_(+over-)')
    predicted_enantiomeric_ratio_plus_over_minus = models.FloatField(db_column='predicted_enantiomeric_ratio_(+over-)')
    seed_address = models.TextField()
    seed_address_alphanumeric = models.TextField()
    seed_plate = models.TextField()
    library_barcode = models.TextField()
    library_ref = models.TextField()
    library_description = models.TextField()
    seed_plate_time = models.TextField()
    pellet_OD = models.FloatField()
    pellet_detected = models.BooleanField()
    pellet_area = models.FloatField()
    pellet_intensity = models.FloatField()
