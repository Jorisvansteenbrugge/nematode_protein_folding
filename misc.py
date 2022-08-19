def strip_pdb_ids(ids):
    return [x.replace('.pdb','') for x in ids]
