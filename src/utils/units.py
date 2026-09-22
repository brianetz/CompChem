HARTREE_TO_KCAL_MOL = 627.509474
HARTREE_TO_EV = 27.211386245988
BOHR_TO_ANGSTROM = 0.529177210903


def hartree_to_kcal(value: float) -> float:
    return value * HARTREE_TO_KCAL_MOL


def delta_hartree_to_kcal(delta: float) -> float:
    return delta * HARTREE_TO_KCAL_MOL


def ev_to_hartree(value: float) -> float:
    return value / HARTREE_TO_EV

