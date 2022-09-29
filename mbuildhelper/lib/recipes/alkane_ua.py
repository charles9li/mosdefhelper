"""mBuild recipe for a generic united-atom alkane chain."""
import mbuild as mb
from mbuildhelper.lib.atoms import CH2sp3UA, CH3UA
from mbuildhelper.lib.molecules import MethaneUA, EthaneUA


class AlkaneUA(mb.Compound):
    """A united-atom alkane which may optionally end with a Port."""

    def __init__(self, n=3, cap_front=True, cap_end=True):
        if n < 1:
            raise ValueError("n must be 1 or more")
        super(AlkaneUA, self).__init__()

        # handle the case of methane and ethane separately
        if n < 3:
            if n == 1:
                if cap_front and cap_end:
                    self.add(MethaneUA(), "chain")
                elif cap_front != cap_end:
                    chain = CH3UA()
                    self.add(chain, "chain")
                    if cap_front:
                        self.add(chain["up"], "down", containment=False)
                    else:
                        self.add(chain["up"], "up", containment=False)
                else:
                    chain = CH2sp3UA()
                    self.add(chain)
                    self.add(chain["down"], "down", containment=False)
                    self.add(chain["up"], "up", containment=False)
            elif n == 2:
                if cap_front and cap_end:
                    self.add(EthaneUA(), "chain")
                elif cap_front != cap_end:
                    chain = CH2sp3UA()
                    self.add(chain, "chain")
                    if cap_front:
                        self.add(CH3UA(), "methyl_front")
                        mb.force_overlap(
                            move_this=self["chain"],
                            from_positions=self["chain"]["up"],
                            to_positions=self["methyl_Front"]["up"]
                        )
                        self.add(chain["down"], "down", containment=False)
                    else:
                        self.add(CH3UA(), "methyl_end")
                        mb.force_overlap(
                            move_this=self["chain"],
                            from_positions=self["methyl_end"]["up"],
                            to_positions=self["chain"]["down"]
                        )
                        self.add(chain["up"], "up", containment=False)
                else:
                    chain = CH2sp3UA()
                    self.add(chain, "chain")
                    self.add(chain["down"], "down", containment=False)
                    self.add(chain["up"], "up", containment=False)

        # handle the general case of n >= 3
        else:
            last_part = None
            if cap_front:
                n -= 1
                this_part = CH3UA()
                self.add(this_part, "monomer[$]")
                last_part = this_part
            if cap_end:
                n -= 1
            for i in range(n):
                this_part = CH2sp3UA()
                self.add(this_part, "monomer[$]")
                if i == 0 and cap_front:
                    mb.force_overlap(
                        move_this=this_part,
                        from_positions=this_part["up"],
                        to_positions=last_part["up"]
                    )
                else:
                    if last_part is not None:
                        mb.force_overlap(
                            move_this=this_part,
                            from_positions=this_part["up"],
                            to_positions=last_part["down"]
                        )
                last_part = this_part
            if cap_end:
                this_part = CH3UA()
                self.add(this_part, "monomer[$]")
                mb.force_overlap(
                    move_this=this_part,
                    from_positions=this_part["up"],
                    to_positions=last_part["down"]
                )
            if not cap_front:
                self.add(self["monomer[0]"]["up"], "up", containment=False)
            if not cap_end:
                self.add(self[f"monomer[{n+int(cap_front)+int(cap_end)-1}]"]["down"], "down", containment=False)

            # end_groups = [None, None]
            # # adjust length of Polymer for absence of methyl terminations
            # if cap_front:
            #     n -= 1
            #     end_groups[0] = CH3UA()
            # if cap_end:
            #     n -= 1
            #     end_groups[1] = CH3UA()

            # chain = mb.recipes.Polymer(monomers=[CH2sp3UA()], end_groups=end_groups)
            # chain.build(n, add_hydrogens=False)
            # self.add(chain, "chain")
            # if not cap_front:
            #     # hoist port to alkane level
            #     self.add(self["chain"]["up"], "up", containment=False)
            # if not cap_end:
            #     # hoist port to alkane level
            #     self.add(self["chain"]["down"], "down", containment=False)


if __name__ == '__main__':
    m = AlkaneUA(n=12)
    m.save("dodecane_ua.mol2", overwrite=True)
