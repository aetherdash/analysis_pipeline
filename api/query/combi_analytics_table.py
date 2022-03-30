from typing import Any, Iterable, List, Optional

import graphene
from promise import Promise
from analysis.models import CombiAnalyticsTable
from federation.directives import key


@key(fields=['id'])
class CombiAnalyticsTableType(graphene.ObjectType):

    class Meta:
        interfaces = (graphene.Node, )
        
    exp_workflow_barcode = graphene.String()
    exp_workflow_name = graphene.String()
    proj_barcode = graphene.String()
    proj_name = graphene.String()
    maldi_run = graphene.String()
    lcms_C18_run = graphene.String()
    lcms_chiral_run = graphene.String()
    maldi_plate_racemic = graphene.String()
    maldi_address_racemic = graphene.String()
    lcms_plate_racemic = graphene.String()
    lcms_address_racemic = graphene.String()
    source_plate_racemic = graphene.String()
    source_address_racemic = graphene.String()
    maldi_plate_plus = graphene.String()
    maldi_address_plus = graphene.String()
    lcms_plate_plus = graphene.String()
    lcms_address_plus = graphene.String()
    source_plate_plus = graphene.String()
    source_address_plus = graphene.String()
    maldi_plate_minus = graphene.String()
    maldi_address_minus = graphene.String()
    lcms_plate_minus = graphene.String()
    lcms_address_minus = graphene.String()
    source_plate_minus = graphene.String()
    source_address_minus = graphene.String()
    ctrl_type = graphene.String()
    exp_condition = graphene.String()
    enzyme_barcode = graphene.String()
    sequence = graphene.String()
    mutations = graphene.String()
    hamming = graphene.Int()
    reference_enzyme = graphene.String()
    substrate_barcode_racemic = graphene.String()
    substrate_barcode_plus = graphene.String()
    substrate_barcode_minus = graphene.String()
    substrate_smiles_racemic = graphene.String()
    substrate_smiles_plus = graphene.String()
    substrate_smiles_minus = graphene.String()
    substrate_concentration_racemic = graphene.Float()
    substrate_concentration_plus = graphene.Float()
    substrate_concentration_minus = graphene.Float()
    predicted_binary_score_racemic = graphene.Float()
    predicted_nonbinary_score_racemic = graphene.Float()
    measured_nonbinary_score_racemic = graphene.Float()
    measured_nonbinary_sum_racemic = graphene.Float()
    measured_conversion_racemic = graphene.Float()
    predicted_binary_score_variance_racemic = graphene.Float()
    predicted_nonbinary_score_variance_racemic = graphene.Float()
    predicted_binary_score_plus = graphene.Float()
    predicted_nonbinary_score_plus = graphene.Float()
    measured_nonbinary_score_plus = graphene.Float()
    measured_conversion_plus = graphene.Float()
    predicted_binary_score_variance_plus = graphene.Float()
    predicted_nonbinary_score_variance_plus = graphene.Float()
    predicted_binary_score_minus = graphene.Float()
    predicted_nonbinary_score_minus = graphene.Float()
    measured_nonbinary_score_minus = graphene.Float()
    measured_conversion_minus = graphene.Float()
    predicted_binary_score_variance_minus = graphene.Float()
    predicted_nonbinary_score_variance_minus = graphene.Float()
    measured_enantiomeric_excess_plus_over_minus = graphene.Float()
    predicted_enantiomeric_excess_plus_over_minus = graphene.Float()
    measured_enantiomeric_ratio_plus_over_minus = graphene.Float()
    predicted_enantiomeric_ratio_plus_over_minus = graphene.Float()
    seed_address = graphene.String()
    seed_address_alphanumeric = graphene.String()
    seed_plate = graphene.String()
    library_barcode = graphene.String()
    library_ref = graphene.String()
    library_description = graphene.String()
    seed_plate_time = graphene.String()
    pellet_OD = graphene.Float()
    pellet_detected = graphene.Boolean()
    pellet_area = graphene.Float()
    pellet_intensity = graphene.Float()

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
    def resolve_maldi_run(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.maldi_run

    @staticmethod
    def resolve_lcms_C18_run(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.lcms_C18_run

    @staticmethod
    def resolve_lcms_chiral_run(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.lcms_chiral_run

    @staticmethod
    def resolve_maldi_plate_racemic(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.maldi_plate_racemic

    @staticmethod
    def resolve_maldi_address_racemic(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.maldi_address_racemic

    @staticmethod
    def resolve_lcms_plate_racemic(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.lcms_plate_racemic

    @staticmethod
    def resolve_lcms_address_racemic(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.lcms_address_racemic

    @staticmethod
    def resolve_source_plate_racemic(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.source_plate_racemic

    @staticmethod
    def resolve_source_address_racemic(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.source_address_racemic

    @staticmethod
    def resolve_maldi_plate_plus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.maldi_plate_plus

    @staticmethod
    def resolve_maldi_address_plus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.maldi_address_plus

    @staticmethod
    def resolve_lcms_plate_plus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.lcms_plate_plus

    @staticmethod
    def resolve_lcms_address_plus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.lcms_address_plus

    @staticmethod
    def resolve_source_plate_plus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.source_plate_plus

    @staticmethod
    def resolve_source_address_plus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.source_address_plus

    @staticmethod
    def resolve_maldi_plate_minus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.maldi_plate_minus

    @staticmethod
    def resolve_maldi_address_minus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.maldi_address_minus

    @staticmethod
    def resolve_lcms_plate_minus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.lcms_plate_minus

    @staticmethod
    def resolve_lcms_address_minus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.lcms_address_minus

    @staticmethod
    def resolve_source_plate_minus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.source_plate_minus

    @staticmethod
    def resolve_source_address_minus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.source_address_minus

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
    def resolve_mutations(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.mutations

    @staticmethod
    def resolve_hamming(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.hamming

    @staticmethod
    def resolve_reference_enzyme(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.reference_enzyme

    @staticmethod
    def resolve_substrate_barcode_racemic(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.substrate_barcode_racemic

    @staticmethod
    def resolve_substrate_barcode_plus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.substrate_barcode_plus

    @staticmethod
    def resolve_substrate_barcode_minus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.substrate_barcode_minus

    @staticmethod
    def resolve_substrate_smiles_racemic(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.substrate_smiles_racemic

    @staticmethod
    def resolve_substrate_smiles_plus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.substrate_smiles_plus

    @staticmethod
    def resolve_substrate_smiles_minus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.substrate_smiles_minus

    @staticmethod
    def resolve_substrate_concentration_racemic(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.substrate_concentration_racemic

    @staticmethod
    def resolve_substrate_concentration_plus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.substrate_concentration_plus

    @staticmethod
    def resolve_substrate_concentration_minus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.substrate_concentration_minus

    @staticmethod
    def resolve_predicted_binary_score_racemic(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.predicted_binary_score_racemic

    @staticmethod
    def resolve_predicted_nonbinary_score_racemic(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.predicted_nonbinary_score_racemic

    @staticmethod
    def resolve_measured_nonbinary_score_racemic(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.measured_nonbinary_score_racemic

    @staticmethod
    def resolve_measured_nonbinary_sum_racemic(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.measured_nonbinary_sum_racemic

    @staticmethod
    def resolve_measured_conversion_racemic(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.measured_conversion_racemic

    @staticmethod
    def resolve_predicted_binary_score_variance_racemic(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.predicted_binary_score_variance_racemic

    @staticmethod
    def resolve_predicted_nonbinary_score_variance_racemic(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.predicted_nonbinary_score_variance_racemic

    @staticmethod
    def resolve_predicted_binary_score_plus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.predicted_binary_score_plus

    @staticmethod
    def resolve_predicted_nonbinary_score_plus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.predicted_nonbinary_score_plus

    @staticmethod
    def resolve_measured_nonbinary_score_plus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.measured_nonbinary_score_plus

    @staticmethod
    def resolve_measured_conversion_plus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.measured_conversion_plus

    @staticmethod
    def resolve_predicted_binary_score_variance_plus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.predicted_binary_score_variance_plus

    @staticmethod
    def resolve_predicted_nonbinary_score_variance_plus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.predicted_nonbinary_score_variance_plus

    @staticmethod
    def resolve_predicted_binary_score_minus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.predicted_binary_score_minus

    @staticmethod
    def resolve_predicted_nonbinary_score_minus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.predicted_nonbinary_score_minus

    @staticmethod
    def resolve_measured_nonbinary_score_minus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.measured_nonbinary_score_minus

    @staticmethod
    def resolve_measured_conversion_minus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.measured_conversion_minus

    @staticmethod
    def resolve_predicted_binary_score_variance_minus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.predicted_binary_score_variance_minus

    @staticmethod
    def resolve_predicted_nonbinary_score_variance_minus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.predicted_nonbinary_score_variance_minus

    @staticmethod
    def resolve_measured_enantiomeric_excess_plus_over_minus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.measured_enantiomeric_excess_plus_over_minus

    @staticmethod
    def resolve_predicted_enantiomeric_excess_plus_over_minus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.predicted_enantiomeric_excess_plus_over_minus

    @staticmethod
    def resolve_measured_enantiomeric_ratio_plus_over_minus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.measured_enantiomeric_ratio_plus_over_minus

    @staticmethod
    def resolve_predicted_enantiomeric_ratio_plus_over_minus(root: None, info: graphene.ResolveInfo, **kwargs) -> Iterable[Any]:
        return root.predicted_enantiomeric_ratio_plus_over_minus

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
