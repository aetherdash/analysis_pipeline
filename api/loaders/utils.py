from collections import defaultdict

from django.apps import apps
from django.conf import settings
from django.db.models import Q

from promise import Promise
from promise.dataloader import DataLoader as BaseDataLoader


def get_bulk_lookup(values, key_fn, many):
    default = list if many else lambda: None
    bulk_lookup = defaultdict(default)
    for value in values:
        key = getattr(value, key_fn)
        if many:
            bulk_lookup[key].append(value)
        else:
            bulk_lookup[key] = value
    return dict(bulk_lookup)


def batch_values(bulk_lookup, batched_keys, many):
    default = [] if many else None
    batched_values = [bulk_lookup.get(key, default) for key in batched_keys]
    return batched_values


def batch_load_model_by_field(app_name, model_name, field_name):
    model_to_load = apps.get_model(app_name, model_name)
    field = model_to_load._meta.get_field(field_name)  # pylint: disable=protected-access
    many = field.many_to_one

    def get_values(batched_keys):
        condition = {f'{field.column}__in':batched_keys}
        return model_to_load.objects.filter(**condition)

    def get_batched_values(batched_keys):
        values = get_values(batched_keys)
        bulk_lookup = get_bulk_lookup(values, field.column, many)
        batched_values = batch_values(bulk_lookup, batched_keys, many)
        return batched_values
    return get_batched_values


def batch_load_model_by_complex_key(app_name, model_name, field_names, filter_operator):
    model_to_load = apps.get_model(app_name, model_name)
    fields = (model_to_load._meta.get_field(field_name) for field_name in field_names)  # pylint: disable=protected-access

    def load(batched_keys):
        queryset = Q()

        for idx, field in enumerate(fields):
            conditions = defaultdict(list)

            for field_tuple in batched_keys:
                conditions[f'{field.column}__in'].append(field_tuple[idx])

            queryset.add( Q(**conditions), filter_operator)

        return model_to_load.objects.filter(queryset)

    def batch_load_fn(batched_keys):
        model_records = load(batched_keys)
        model_map = {}
        for model_record in model_records:
            field_values = tuple( getattr(model_record, name) for name in field_names )
            model_map[ field_values ] = model_record

        return [model_map.get(key) for key in batched_keys]

    return batch_load_fn


def batch_load_fancy(app_name, model_name):
    model_to_load = apps.get_model(app_name, model_name)
    def load(primary_keys):
        return model_to_load.objects.filter(object_id__in=primary_keys)

    def batch_load_fn(model_ids):
        model_records = load(model_ids)
        model_map = {}

        for model_record in model_records:
            model_map[model_record.object_id] = model_record

        return [model_map.get(model_id) for model_id in model_ids]
    return batch_load_fn


def batch_load_primary_key(app_name, model_name):
    model_to_load = apps.get_model(app_name, model_name)
    def load(primary_keys):
        return model_to_load.objects.filter(id__in=primary_keys)

    def batch_load_fn(model_ids):
        model_records = load(model_ids)
        model_map = {}
        for model_record in model_records:
            model_map[model_record.id] = model_record
        return [model_map.get(model_id) for model_id in model_ids]
    return batch_load_fn


def batch_load_foreign_key(app_name, model_name, foreign_key_name):
    model_to_load = apps.get_model(app_name, model_name)
    # key_field = model_to_load._meta.get_field(foreign_key_name)
    def load(foreign_keys):
        condition = {f'{foreign_key_name}_id__in':foreign_keys}
        return model_to_load.objects.filter(**condition)

    def batch_load_fn(model_ids):
        model_records = load(model_ids)
        model_map = defaultdict(list)
        for model in model_records:
            foreign_key_id = getattr(model, f'{foreign_key_name}_id')
            model_map[foreign_key_id].append(model)
        return [model_map[model_id] for model_id in model_ids]
    return batch_load_fn


def batch_load_many_to_many_key(app_name, model_name, foreign_key_name, desired_key_name):
    model_to_load = apps.get_model(app_name, model_name)
    # key_field = model_to_load._meta.get_field(foreign_key_name)
    def load(foreign_keys):
        condition = {f'{foreign_key_name}_id__in':foreign_keys}
        return model_to_load.objects.filter(**condition)

    def batch_load_fn(model_ids):
        model_records = load(model_ids)
        model_map = defaultdict(list)
        for model in model_records:
            foreign_key_id = getattr(model, f'{foreign_key_name}_id')
            model_record = getattr(model, f'{desired_key_name}')
            model_map[foreign_key_id].append(model_record)
        return [model_map[model_id] for model_id in model_ids]
    return batch_load_fn


class DataLoader(BaseDataLoader):
    def __init__(self, *args, **kwargs):
        load_fn = kwargs.pop('batch_load_fn', None)
        if not load_fn:
            load_fn = args[0]
            args = (None,) + args[1:]
        self.load_fn = load_fn
        self.max_batch_size = settings.LOADER_BATCH_SIZE

        super().__init__(*args, **kwargs)

    def batch_load_fn(self, keys): # pylint: disable=method-hidden
        return Promise.resolve(self.load_fn(keys))
