"""An ethane molecule for united-atom force fields."""
import mbuild as mb
from mbuildhelper.lib.atoms import CH3UA


class EthaneUA(mb.Compound):
    """An ethane molecule for united-atom force fields.

    Connect two united-atom methyl groups to form an ethane.
    """

    def __init__(self):
        super(EthaneUA, self).__init__()

        self.add(CH3UA(), "methyl1")
        self.add(CH3UA(), "methyl2")
        mb.force_overlap(
            self["methyl1"], self["methyl1"]["up"], self["methyl2"]["up"]
        )


if __name__ == '__main__':
    m = EthaneUA()
    m.save("ethane_ua.mol2", overwrite=True)
