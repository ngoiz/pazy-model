from structure import PazyStructure
import sharpy.utils.settings as settings

class PazyWing:

    model_settings_default = dict()
    model_settings_types = dict()
    model_settings_options = dict()
    model_settings_description = dict()

    model_settings_default['skin_on'] = False
    model_settings_types['skin_on'] = 'bool'

    model_settings_default['discretisation_method'] = 'michigan'
    model_settings_types['discretisation_method'] = 'str'
    model_settings_options['discretisation_method'] = ['michigan', 'even', 'fine_root_tip']

    model_settings_types['num_elem'] = 'int'
    model_settings_default['num_elem'] = 2
    model_settings_description['num_elem'] = 'If discretisation is ``michigan`` then it corresponds to how many ' \
                                             'times ' \
                                             'the original number of elements is replicated, else is the number of ' \
                                             'elements'

    def __init__(self, case_name, case_route='./', in_settings=None):
        self.case_name = case_name
        self.case_route = case_route

        if in_settings is not None:
            self.settings = in_settings
        else:
            self.settings = dict()

        settings.to_custom_types(self.settings, self.model_settings_types, self.model_settings_default,
                                 self.model_settings_options, no_ctype=True)

        self.structure = None

    def generate_structure(self):
        pazy = PazyStructure(**self.settings)
        pazy.generate()

        self.structure = pazy

    def save_files(self):
        self.structure.save_fem_file(case_name=self.case_name, case_route=self.case_route)

# load michigan's results

if __name__ == '__main__':
    import sharpy.sharpy_main
    # for n_elem in [16, 32, 64, 92, 128]:
    #     pazy = PazyWing()
    #     case_name = 'modal_inertia_even_n{}'.format(n_elem)
    #     output_folder = './output/'
    #     pazy.generate_structure(n_elem, case_name=case_name)
    #     pazy.structure.create_modal_simulation(case_name=case_name, output_folder=output_folder)
    #     sharpy.sharpy_main.main(['', case_name + '.sharpy'])

    n_elem = 32
    # pazy = PazyWing()
    case_name = 'modal_inertia_even_n{}'.format(n_elem)
    output_folder = './output/'
    # pazy.generate_structure(n_elem, case_name=case_name)
    # pazy.structure.create_modal_simulation(case_name=case_name, output_folder=output_folder)
    # sharpy.sharpy_main.main(['', case_name + '.sharpy'])

    model_settings = {'skin_on': 'off',
                      'discretisation_method': 'michigan',
                      'num_elem': 2}

    pazy = PazyWing(case_name, model_settings)
    pazy.generate_structure()
    # pazy.structure.create_modal_simulation(case_name=case_name, output_folder=output_folder)

    pazy.structure.add_lumped_mass((10, 0))
    pazy.save_files()
    # sharpy.sharpy_main.main(['', case_name + '.sharpy'])

