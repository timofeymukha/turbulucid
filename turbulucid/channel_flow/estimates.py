__all__ = ["dean_reb_from_retau", "dean_retau_from_reb", "dean_rec_from_reb",
           "dean_reb_from_rec", "dean_rec_from_retau", "dean_retau_from_rec"]


def dean_retau_from_reb(reB):
    return 0.175*reB**0.875


def dean_rec_from_reb(reB):
    return 1.27*reB**0.988


def dean_reb_from_retau(reTau):
    return 1/0.175*reTau**(1/0.875)


def dean_reb_from_rec(reC):
    return 1/1.27*reC**(1/0.988)


def dean_retau_from_rec(reC):
    reB = dean_reb_from_rec(reC)
    return dean_retau_from_reb(reB)


def dean_rec_from_retau(reTau):
    reB = dean_reb_from_retau(reTau)
    return dean_rec_from_reb(reB)

