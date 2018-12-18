def run_cmd(cmd):
    """
    Run given CMD
    :param cmd: CMD (str)
    :return: cmd (str), status (int), stdout + stderr (str)
    """
    import subprocess
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    p_stdout = p.stdout.read().decode()
    p_comm   = p.communicate()[0]
    p_status = p.returncode
    return cmd, p_status, p_stdout

def proc_taxon_query_results(fpath):
    import re, pandas
    header  = ["taxon_id", "taxon_name", "taxon_rank", "lineage"]
    ranks   = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]
    for rank in ranks:
        header += ['taxon_%s_id' % rank, 'taxon_%s_name' % rank]
    taxa = []
    with open(fpath, 'r') as ifile:
        for line in ifile:
            info = dict.fromkeys(header,'NA')
            line = line.rstrip('\n')
            for field in line.split('\t'):
                if re.search('\|',field):
                    tid, tname, trank = field.split('|')
                    if trank == 'no rank':
                        trank = 'NA'
                    elif "taxon_%s_name" % trank in header:
                        info["taxon_%s_name" % trank] = tname
                        info["taxon_%s_id" % trank] = tid
                    if info['taxon_id'] == 'NA':
                        info['taxon_id'] = tid
                        info['taxon_name'] = tname
                        info['taxon_rank'] = trank
                else: # lineage
                    info['lineage'] = field
            assert info['taxon_id'] != "NA" and info['taxon_name'] != "NA"
            taxa.append(dict([(k, info[k]) for k in header]))
    taxa = pandas.DataFrame(taxa)
    return taxa

def aggr_taxa_for_xml(node, df, rank, ranks, labels):
    import xml.etree.cElementTree as ET
    # print('Rank %s, %d entries' % (rank, df.shape[0]))
    for taxon, taxon_df in df.groupby(by=rank):
        taxon_node = ET.SubElement(node, 'node')
        taxon_node.attrib['name'] = taxon
        taxon_count = ET.SubElement(taxon_node, 'count')
        for db in labels:
            if db in set(taxon_df.index.get_level_values(0)):
                db_count = ET.SubElement(taxon_count, 'val').text = str(taxon_df.loc[db,:].shape[0])
            else:
                db_count = ET.SubElement(taxon_count, 'val').text = '0'
        if ranks.index(rank) < len(ranks) - 1:
            taxon_node = aggr_taxa_for_xml(taxon_node, taxon_df, ranks[ranks.index(rank) + 1], ranks, labels)
    return node
