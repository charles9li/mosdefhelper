"""United-atom acrylate and methacrylate moieties."""
__all__ = ['AcrylateUA', 'MethacrylateUA']

import mbuild as mb

from mbuildhelper.lib.atoms import Csp2, Csp3, CHsp2UA, CHsp3UA, CH2sp2UA, CH2sp3UA, CH3UA
from mbuildhelper.lib.moieties import Ester


class AcrylateUA(mb.Compound):

    def __init__(self, cap_front=True, cap_end=True, methyl=False, **kwargs):
        super(AcrylateUA, self).__init__(**kwargs)

        # first handle case of single unpolymerized (meth)acrylate molecule
        if cap_front and cap_end:
            first_carbon = CH2sp2UA()
            self.add(first_carbon)
            first_carbon_port = first_carbon['up']
            if methyl:
                second_carbon = Csp2()
            else:
                second_carbon = CHsp2UA()
            self.add(second_carbon)
            mb.force_overlap(
                move_this=second_carbon,
                from_positions=second_carbon['up'],
                to_positions=first_carbon_port
            )
            second_carbon_port = second_carbon['down']
            if methyl:
                methyl_carbon = CH3UA()
                self.add(methyl_carbon)
                mb.force_overlap(
                    move_this=methyl_carbon,
                    from_positions=methyl_carbon['up'],
                    to_positions=second_carbon['left']
                )
        # all other cases
        else:
            # add first carbon
            if cap_front:
                first_carbon = CH3UA()
                self.add(first_carbon)
                first_carbon_port = first_carbon['up']
            else:
                first_carbon = CH2sp3UA()
                self.add(first_carbon)
                first_carbon_port = first_carbon["down"]
                # hoist port to acrylate level
                self.add(first_carbon["up"], "up", containment=False)

            # add second carbon
            if cap_end:
                if methyl:
                    second_carbon = CHsp3UA()
                    second_carbon_port = second_carbon['left']
                    methyl_port = second_carbon['down']
                else:
                    second_carbon = CH2sp3UA()
                    second_carbon_port = second_carbon['down']
                self.add(second_carbon)
            else:
                if methyl:
                    second_carbon = Csp3()
                    methyl_port = second_carbon['right']
                else:
                    second_carbon = CHsp3UA()
                self.add(second_carbon)
                second_carbon_port = second_carbon["left"]
                # hoist port to acrylate level
                self.add(second_carbon["down"], "down", containment=False)

            # bond first carbon to second carbon
            mb.force_overlap(
                move_this=second_carbon,
                from_positions=second_carbon["up"],
                to_positions=first_carbon_port
            )

            # bond methyl
            if methyl:
                methyl_carbon = CH3UA()
                self.add(methyl_carbon)
                mb.force_overlap(
                    move_this=methyl_carbon,
                    from_positions=methyl_carbon['up'],
                    to_positions=methyl_port
                )

        # add ester group
        self.add(Ester(), label="ester")
        mb.force_overlap(
            move_this=self['ester'],
            from_positions=self['ester']['R1'],
            to_positions=second_carbon_port
        )
        self.add(self['ester']['R2'], 'R', containment=False)   # hoist port to (meth)acrylate level


class MethacrylateUA(AcrylateUA):

    def __init__(self, cap_front=True, cap_end=True, **kwargs):
        super(MethacrylateUA, self).__init__(cap_front=cap_front, cap_end=cap_end, methyl=True, **kwargs)


if __name__ == '__main__':
    m = AcrylateUA(cap_front=True, cap_end=True)
    m.save("acrylate_ua.mol2", overwrite=True)
