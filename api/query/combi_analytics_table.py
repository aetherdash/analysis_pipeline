from typing import Any, Iterable, List, Optional

import graphene
from promise import Promise
from analysis.models import CombiAnalyticsTable
from federation.directives import key


@key(fields=['id'])
class CombiAnalyticsTableType(graphene.ObjectType):

    class Meta:
        interfaces = (graphene.Node, )
    
    address = graphene.String()
    run = graphene.String()
    plate = graphene.String()
    dev_or_prod = graphene.String()
    exp_workflow_barcode = graphene.String()
    exp_workflow_name = graphene.String()
    proj_barcode = graphene.String()
    proj_name = graphene.String()
    ctrl_type = graphene.String()
    exp_condition = graphene.String()
    enzyme_barcode = graphene.String()
    sequence = graphene.String()
    hamming = graphene.Int()
    mutations = graphene.String()
    reference_enzyme = graphene.String()
    enzyme_concentration = graphene.Float()
    enzyme_unit = graphene.String()
    enzyme_class = graphene.String()
    sequence_qc = graphene.String()
    substrate_barcode = graphene.String()
    substrate_concentration = graphene.Float()
    substrate_unit = graphene.String()
    substrate_smiles = graphene.String()
    substrate_mz = graphene.String()
    product_smiles = graphene.String()
    product_mz = graphene.Float()
    pk_0_mz = graphene.Float()
    isotope_score = graphene.Float()
    pk_sub_raw = graphene.Float()
    pk_sub_raw_ref = graphene.Float()
    pk_prod_raw = graphene.Float()
    pk_prod_raw_ref = graphene.Float()
    pk_0_raw = graphene.Float()
    pk_0_raw_ref = graphene.Float()
    pk_mtx_172_raw = graphene.Float()
    pk_mtx_190_raw = graphene.Float()
    pk_mtx_379_raw = graphene.Float()
    pk_sub_stdz_379 = graphene.Float()
    pk_sub_stdz_190 = graphene.Float()
    pk_sub_stdz_172 = graphene.Float()
    pk_prod_stdz_379 = graphene.Float()
    pk_prod_stdz_190 = graphene.Float()
    pk_prod_stdz_172 = graphene.Float()
    pk_0_stdz_379 = graphene.Float()
    pk_0_stdz_190 = graphene.Float()
    pk_0_stdz_172 = graphene.Float()
    ptr = graphene.Float()
    psr = graphene.Float()
    tstat_pk_sub_stdz_379 = graphene.Float()
    tstat_pk_sub_stdz_190 = graphene.Float()
    tstat_pk_sub_stdz_172 = graphene.Float()
    tstat_pk_prod_stdz_379 = graphene.Float()
    tstat_pk_prod_stdz_190 = graphene.Float()
    tstat_pk_prod_stdz_172 = graphene.Float()
    tstat_pk_0_stdz_379 = graphene.Float()
    tstat_pk_0_stdz_190 = graphene.Float()
    tstat_pk_0_stdz_172 = graphene.Float()
    tstat_ptr = graphene.Float()
    tstat_psr = graphene.Float()
    tstat_multi = graphene.Float()
    spectra_ids_ref = graphene.String()
    spectra_ids = graphene.String()
    succinimide_bool = graphene.Boolean()
    tertiary_amide_bool = graphene.Boolean()
    repo_sha = graphene.String()
    parameters = graphene.String()
    rxn_smarts = graphene.String()
    retro_rules_id = graphene.String()
    true_score_binary = graphene.Float()
    true_score_nonbinary = graphene.Float()
    model_score_binary = graphene.Float()
    model_score_nonbinary = graphene.Float()
    timestamp = graphene.String()
    source_plate = graphene.String()
    source_address = graphene.String()
    expression_condition = graphene.String()
    screen_condition = graphene.String()
    acq_condition = graphene.String()
    seed_address = graphene.String()
    seed_address_alphanumeric = graphene.String()
    seed_plate = graphene.String()
    main_plate = graphene.String()
    rxn_plate = graphene.String()
    library_barcode = graphene.String()
    library_ref = graphene.String()
    library_description = graphene.String()
    seed_plate_time = graphene.String()
    main_plate_time = graphene.String()
    rxn_plate_time = graphene.String()
    plate_time = graphene.String()
    pellet_OD = graphene.Float()
    pellet_detected = graphene.Boolean()
    pellet_area = graphene.Float()
    pellet_intensity = graphene.Float()
 
    @staticmethod
    def resolve_address(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.address
        
    @staticmethod
    def resolve_run(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.run
        
    @staticmethod
    def resolve_plate(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.plate
     
    @staticmethod
    def resolve_dev_or_prod(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.dev_or_prod
     
    @staticmethod
    def resolve_exp_workflow_barcode(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.exp_workflow_barcode
     
    @staticmethod
    def resolve_exp_workflow_name(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.exp_workflow_name
     
    @staticmethod
    def resolve_proj_barcode(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.proj_barcode
     
    @staticmethod
    def resolve_proj_name(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.proj_name
     
    @staticmethod
    def resolve_ctrl_type(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.ctrl_type
     
    @staticmethod
    def resolve_exp_condition(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.exp_condition
     
    @staticmethod
    def resolve_enzyme_barcode(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.enzyme_barcode
     
    @staticmethod
    def resolve_sequence(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.sequence

    @staticmethod
    def resolve_hamming(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.hamming
     
    @staticmethod
    def resolve_mutations(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.mutations
     
    @staticmethod
    def resolve_reference_enzyme(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.reference_enzyme

    @staticmethod
    def resolve_enzyme_concentration(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.enzyme_concentration
     
    @staticmethod
    def resolve_enzyme_unit(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.enzyme_unit
     
    @staticmethod
    def resolve_enzyme_class(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.enzyme_class
     
    @staticmethod
    def resolve_sequence_qc(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.sequence_qc
     
    @staticmethod
    def resolve_substrate_barcode(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.substrate_barcode
    
    @staticmethod
    def resolve_substrate_concentration(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.substrate_concentration
     
    @staticmethod
    def resolve_substrate_unit(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.substrate_unit
     
    @staticmethod
    def resolve_substrate_smiles(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.substrate_smiles
     
    @staticmethod
    def resolve_substrate_mz(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.substrate_mz
     
    @staticmethod
    def resolve_product_smiles(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.product_smiles
        
    @staticmethod
    def resolve_product_mz(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.product_mz
     
    @staticmethod
    def resolve_pk_0_mz(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pk_0_mz
     
    @staticmethod
    def resolve_isotope_score(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.isotope_score
     
    @staticmethod
    def resolve_pk_sub_raw(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pk_sub_raw
     
    @staticmethod
    def resolve_pk_sub_raw_ref(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pk_sub_raw_ref
     
    @staticmethod
    def resolve_pk_prod_raw(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pk_prod_raw
     
    @staticmethod
    def resolve_pk_prod_raw_ref(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pk_prod_raw_ref
     
    @staticmethod
    def resolve_pk_0_raw(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pk_0_raw
     
    @staticmethod
    def resolve_pk_0_raw_ref(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pk_0_raw_ref
     
    @staticmethod
    def resolve_pk_mtx_172_raw(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pk_mtx_172_raw
     
    @staticmethod
    def resolve_pk_mtx_190_raw(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pk_mtx_190_raw
     
    @staticmethod
    def resolve_pk_mtx_379_raw(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pk_mtx_379_raw
     
    @staticmethod
    def resolve_pk_sub_stdz_379(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pk_sub_stdz_379
     
    @staticmethod
    def resolve_pk_sub_stdz_190(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pk_sub_stdz_190
     
    @staticmethod
    def resolve_pk_sub_stdz_172(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pk_sub_stdz_172
     
    @staticmethod
    def resolve_pk_prod_stdz_379(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pk_prod_stdz_379
     
    @staticmethod
    def resolve_pk_prod_stdz_190(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pk_prod_stdz_190
     
    @staticmethod
    def resolve_pk_prod_stdz_172(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pk_prod_stdz_172
     
    @staticmethod
    def resolve_pk_0_stdz_379(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pk_0_stdz_379
     
    @staticmethod
    def resolve_pk_0_stdz_190(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pk_0_stdz_190
     
    @staticmethod
    def resolve_pk_0_stdz_172(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pk_0_stdz_172
     
    @staticmethod
    def resolve_ptr(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.ptr
     
    @staticmethod
    def resolve_psr(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.psr
     
    @staticmethod
    def resolve_tstat_pk_sub_stdz_379(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.tstat_pk_sub_stdz_379
     
    @staticmethod
    def resolve_tstat_pk_sub_stdz_190(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.tstat_pk_sub_stdz_190
     
    @staticmethod
    def resolve_tstat_pk_sub_stdz_172(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.tstat_pk_sub_stdz_172
     
    @staticmethod
    def resolve_tstat_pk_prod_stdz_379(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.tstat_pk_prod_stdz_379
     
    @staticmethod
    def resolve_tstat_pk_prod_stdz_190(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.tstat_pk_prod_stdz_190
     
    @staticmethod
    def resolve_tstat_pk_prod_stdz_172(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.tstat_pk_prod_stdz_172
     
    @staticmethod
    def resolve_tstat_pk_0_stdz_379(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.tstat_pk_0_stdz_379
     
    @staticmethod
    def resolve_tstat_pk_0_stdz_190(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.tstat_pk_0_stdz_190
     
    @staticmethod
    def resolve_tstat_pk_0_stdz_172(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.tstat_pk_0_stdz_172
     
    @staticmethod
    def resolve_tstat_ptr(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.tstat_ptr
     
    @staticmethod
    def resolve_tstat_psr(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.tstat_psr
     
    @staticmethod
    def resolve_tstat_multi(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.tstat_multi
     
    @staticmethod
    def resolve_spectra_ids_ref(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.spectra_ids_ref
     
    @staticmethod
    def resolve_spectra_ids(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.spectra_ids

    @staticmethod
    def resolve_succinimide_bool(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.succinimide_bool
     
    @staticmethod
    def resolve_tertiary_amide_bool(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.tertiary_amide_bool
     
    @staticmethod
    def resolve_repo_sha(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.repo_sha
     
    @staticmethod
    def resolve_parameters(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.parameters
     
    @staticmethod
    def resolve_rxn_smarts(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.rxn_smarts
     
    @staticmethod
    def resolve_retro_rules_id(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.retro_rules_id

    @staticmethod
    def resolve_true_score_binary(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.true_score_binary
     
    @staticmethod
    def resolve_true_score_nonbinary(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.true_score_nonbinary
     
    @staticmethod
    def resolve_model_score_binary(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.model_score_binary
     
    @staticmethod
    def resolve_model_score_nonbinary(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.model_score_nonbinary

    @staticmethod
    def resolve_timestamp(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.timestamp
     
    @staticmethod
    def resolve_source_plate(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.source_plate
     
    @staticmethod
    def resolve_source_address(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.source_address
     
    @staticmethod
    def resolve_expression_condition(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.expression_condition
     
    @staticmethod
    def resolve_screen_condition(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.screen_condition
     
    @staticmethod
    def resolve_acq_condition(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.acq_condition
     
    @staticmethod
    def resolve_seed_address(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.seed_address
     
    @staticmethod
    def resolve_seed_address_alphanumeric(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.seed_address_alphanumeric
     
    @staticmethod
    def resolve_seed_plate(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.seed_plate
     
    @staticmethod
    def resolve_main_plate(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.main_plate
     
    @staticmethod
    def resolve_rxn_plate(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.rxn_plate
     
    @staticmethod
    def resolve_library_barcode(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.library_barcode
     
    @staticmethod
    def resolve_library_ref(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.library_ref
     
    @staticmethod
    def resolve_library_description(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.library_description
     
    @staticmethod
    def resolve_seed_plate_time(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.seed_plate_time
     
    @staticmethod
    def resolve_main_plate_time(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.main_plate_time
     
    @staticmethod
    def resolve_rxn_plate_time(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.rxn_plate_time
     
    @staticmethod
    def resolve_plate_time(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.plate_time

    @staticmethod
    def resolve_pellet_OD(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pellet_OD
     
    @staticmethod
    def resolve_pellet_detected(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pellet_detected
     
    @staticmethod
    def resolve_pellet_area(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pellet_area
     
    @staticmethod
    def resolve_pellet_intensity(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.pellet_intensity

    @classmethod
    def is_type_of(cls, root: Any, info: graphene.ResolveInfo) -> bool:
        return isinstance(root, CombiAnalyticsTable)


    @classmethod
    def get_node(cls, info: graphene.ResolveInfo, decoded_id: str) -> Promise[Optional[Any]]:
        key = int(decoded_id)
        return info.context.loaders.analysis.load(key)


class CombiAnalyticsTableConnection(graphene.Connection):

    class Meta:
        node = CombiAnalyticsTableType


class CombiAnalyticsTableQuery(graphene.ObjectType):

    analysis = graphene.Node.Field(CombiAnalyticsTableType)
    all_analysis = graphene.ConnectionField(CombiAnalyticsTableConnection)

    @staticmethod
    def resolve_all_analysis(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return CombiAnalyticsTable.objects.all()
