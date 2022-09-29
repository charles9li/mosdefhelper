"""A CH4 (methane) pseudoatom for united-atom force fields."""
import mbuild as mb


class CH4UA(mb.Compound):
    """A united-atom methane."""

    def __init__(self):
        super(CH4UA, self).__init__()
        self.add(mb.Particle(name='_CH4', pos=[0, 0, 0]))


if __name__ == '__main__':
    m = CH4UA()
    m.save("ch4_ua.mol2", overwrite=True)
