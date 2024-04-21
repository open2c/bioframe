import pandas as pd


def parse_gtf_attributes(attrs, kv_sep="=", item_sep=";", quotechar='"', **kwargs):
    item_lists = attrs.str.split(item_sep)
    item_lists = item_lists.apply(
        lambda items: [item.strip().split(kv_sep) for item in items]
    )
    stripchars = quotechar + " "
    item_lists = item_lists.apply(
        lambda items: [
            [x.strip(stripchars) for x in item] for item in items if len(item) == 2
        ]
    )
    kv_records = item_lists.apply(dict)
    return pd.DataFrame.from_records(kv_records, **kwargs)
