SO_LABELS = ('X1SPNORB', 'Y1SPNORB', 'Z1SPNORB')

def get_ls(cmo, symorb):
    orbitals = get_orbitals(cmo, symorb)
    spin_orbit_matrices = (prop.read(m) for m in SO_LABELS)
    return [p.T*ls*p for p in spin_orbit_matrices]

def get_orbitals(cmo, symorb):
    indices = get_orbital_indices(cmo, symorb)
    return cmo.unblock()[:, (9, 15, 23)]


def get_orbital_indices(cmo, symorb):
    indices = (sum(cmo.nrow[:sym-1]) for sym in symorb)
    return tuple(indices)

    
