import unittest

import mbuild as mb

from mbuildhelper.lib.moieties import AcrylateUA, MethacrylateUA
from mbuildhelper.lib.recipes import AlkaneUA


class TestAcrylateUA(unittest.TestCase):

    def test_acrylate_ua(self):
        # acrylate with both ends capped
        acrylate_ua = AcrylateUA()
        acrylate_ua.save("acrylate_ua.mol2", overwrite=True)
        self.assertEqual(1, len(acrylate_ua.available_ports()))
        self.assertSetEqual({"['ester']['R2']", "['R']"}, acrylate_ua.all_ports()[0].access_labels)

        # acrylate with front capped
        acrylate_ua = AcrylateUA(cap_front=True, cap_end=False)
        acrylate_ua.save("acrylate_ua_front_cap.mol2", overwrite=True)
        self.assertEqual(2, len(acrylate_ua.available_ports()))

        # acrylate with end capped
        acrylate_ua = AcrylateUA(cap_front=False, cap_end=True)
        acrylate_ua.save("acrylate_ua_end_cap.mol2", overwrite=True)
        self.assertEqual(2, len(acrylate_ua.available_ports()))

        # acrylate with both sides uncapped
        acrylate_ua = AcrylateUA(cap_front=False, cap_end=False)
        acrylate_ua.save("acrylate_ua_no_cap.mol2", overwrite=True)
        self.assertEqual(3, len(acrylate_ua.available_ports()))
        
    def test_methacrylate_ua(self):
        # methacrylate with both ends capped
        methacrylate_ua = MethacrylateUA()
        methacrylate_ua.save("methacrylate_ua.mol2", overwrite=True)
        self.assertEqual(1, len(methacrylate_ua.available_ports()))
        self.assertSetEqual({"['ester']['R2']", "['R']"}, methacrylate_ua.all_ports()[0].access_labels)

        # methacrylate with front capped
        methacrylate_ua = MethacrylateUA(cap_front=True, cap_end=False)
        methacrylate_ua.save("methacrylate_ua_front_cap.mol2", overwrite=True)
        self.assertEqual(2, len(methacrylate_ua.available_ports()))

        # # acrylate with end capped
        methacrylate_ua = MethacrylateUA(cap_front=False, cap_end=True)
        methacrylate_ua.save("methacrylate_ua_end_cap.mol2", overwrite=True)
        self.assertEqual(2, len(methacrylate_ua.available_ports()))

        # acrylate with both sides uncapped
        methacrylate_ua = MethacrylateUA(cap_front=False, cap_end=False)
        methacrylate_ua.save("methacrylate_ua_no_cap.mol2", overwrite=True)
        self.assertEqual(3, len(methacrylate_ua.available_ports()))

    def test_butyl_acrylate(self):
        butyl_acrylate = mb.Compound()
        butyl_acrylate.add(AcrylateUA(cap_front=False, cap_end=False), label="acrylate")
        butyl_acrylate.add(AlkaneUA(n=4, cap_front=False, cap_end=True), label="butyl")
        mb.force_overlap(
            move_this=butyl_acrylate['butyl'],
            from_positions=butyl_acrylate['butyl']['up'],
            to_positions=butyl_acrylate['acrylate']['R']
        )
        butyl_acrylate.save("butyl_acrylate_ua.mol2", overwrite=True)


if __name__ == '__main__':
    unittest.main()
