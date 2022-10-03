import os
from pkg_resources import resource_filename

import foyer


class ForceField(foyer.Forcefield):

    def __init__(self, forcefield_files=None, name=None, validation=True, debug=False):
        try:
            super(ForceField, self).__init__(forcefield_files=forcefield_files, name=name, validation=validation, debug=debug)
        except (FileNotFoundError, OSError) as err:
            ff_filepath = resource_filename("foyerhelper", "forcefields")
            if forcefield_files is not None:
                if isinstance(forcefield_files, str):
                    forcefield_files = [forcefield_files]
                forcefield_files = [os.path.join(ff_filepath, "xml", fff) for fff in forcefield_files]
                super(ForceField, self).__init__(forcefield_files=forcefield_files)
            elif name is not None:
                super(ForceField, self).__init__(forcefield_files=os.path.join(ff_filepath, "xml", name))
            else:
                raise err

