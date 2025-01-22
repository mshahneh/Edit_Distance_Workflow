def get_edges_index(mol, match):
    """for a mol and a match, returns what edges are in the match
    
    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        the molecule
    match : list
        list of index of the atoms that are in the match
    
    """
    binary_edge = 0
    for i in range(len(match)):
        for j in range(i+1, len(match)):
            bond = mol.GetBondBetweenAtoms(match[i], match[j])
            if bond is not None:
                binary_edge += 1<<bond.GetIdx()
    return binary_edge


def get_sections(dp, par, final):
    item = final
    atom_group = dict()
    index = 0
    while par[item][0] is not None:
        for atom in par[item][0][0]:
            atom_group[atom] = index
        index += 1
        item = par[item][1]
    return atom_group
    

def bfs(all_motif_matches, final, top_limit=None):
    dp = dict()
    par = dict()
    forward_dp = dict()
    forward_par = dict()
    
    if final in all_motif_matches:
        return 1, {atom: 0 for atom in all_motif_matches[final][0]}
    
    for motif_cover in all_motif_matches:
        move_back = final - motif_cover
        forward_dp[move_back] = 1
        forward_par[move_back] = (all_motif_matches[motif_cover], motif_cover)
        
    dp[0] = 0
    par[0] = (None, None)

    just_updated = set(dp.keys())
    count = 1
    while True:
        new_just_updated = set()
        for edge_status in just_updated:
            for motif_cover in all_motif_matches:
                if edge_status | motif_cover == edge_status + motif_cover:
                    new_edge_status = edge_status + motif_cover
                    if new_edge_status not in dp:
                        dp[new_edge_status] = dp[edge_status] + 1
                        par[new_edge_status] = (all_motif_matches[motif_cover], edge_status)
                        new_just_updated.add(new_edge_status)
                    elif dp[new_edge_status] > dp[edge_status] + 1:
                        dp[new_edge_status] = dp[edge_status] + 1
                        par[new_edge_status] = (all_motif_matches[motif_cover], edge_status)
                        new_just_updated.add(new_edge_status)
                
                    if new_edge_status == final:
                        break
                    
                    if new_edge_status in forward_dp:
                        dp[final] = dp[new_edge_status] + forward_dp[new_edge_status]
                        while new_edge_status != final:
                            next_motif, next_edge = forward_par[new_edge_status]
                            next_edge = new_edge_status + next_edge
                            dp[next_edge] = dp[new_edge_status] + 1
                            par[next_edge] = (next_motif, new_edge_status)
                            new_edge_status = next_edge
                        break
                        
            
            if final in dp:
                break
            
        if final in dp:
            break
        just_updated = new_just_updated.copy()
        if len(just_updated) == 0:
            break
        
        count += 1
        if top_limit is not None and count >= top_limit:
            break
    
    if final not in dp:
        return top_limit + 1, None
    else:
        return dp[final], get_sections(dp, par, final)
    

def get_min_motif_to_cover(mol, motif_vocab, top_limit=None):
    # first, create a dictionary of all the possible matches
    all_motif_matches = dict()
    for motif in motif_vocab:
        matches = mol.GetSubstructMatches(motif)
        for match in matches:
            match_index = get_edges_index(mol, match)
            if match_index not in all_motif_matches:
                all_motif_matches[match_index] = (match, motif)
    
    # now, we need to find the minimum set of motifs that cover all the edges
    final = 2**mol.GetNumBonds() - 1
    res, sections = bfs(all_motif_matches, final, top_limit)
    return res, sections
    