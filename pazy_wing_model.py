from structure import PazyStructure
import sharpy.utils.settings as settings
import aero
import configobj
import sharpy.utils.cout_utils as cout


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

    model_settings_types['surface_m'] = 'int'
    model_settings_default['surface_m'] = 4

    model_settings_types['num_surfaces'] = 'int'
    model_settings_default['num_surfaces'] = 2

    def __init__(self, case_name, case_route='./', in_settings=None):
        self.case_name = case_name
        self.case_route = case_route

        cout.start_writer()

        if in_settings is not None:
            self.settings = in_settings
        else:
            self.settings = dict()

        settings.to_custom_types(self.settings, self.model_settings_types, self.model_settings_default,
                                 self.model_settings_options, no_ctype=True)

        self.structure = None
        self.aero = None

        self.config = configobj.ConfigObj()
        self.config.filename = self.case_route + '/' + self.case_name + '.sharpy'

    def save_config(self):
        self.config.write()

    def _get_ea_reference_line(self):

        if self.settings['skin_on']:
            ea = 0.5310
        else:
            ea = 0.4410

        return ea

    def generate_structure(self):
        pazy = PazyStructure(**self.settings)
        pazy.generate()

        self.structure = pazy

    def generate_aero(self):
        pazy_aero = aero.PazyAero(main_ea=self._get_ea_reference_line(),
                                  pazy_structure=self.structure,
                                  **self.settings)

        pazy_aero.generate_aero()

        self.aero = pazy_aero

    def create_aeroelastic_model(self):
        self.generate_structure()
        self.structure.mirror_wing()
        self.generate_aero()

    def save_files(self):
        self.structure.save_fem_file(case_name=self.case_name, case_route=self.case_route)
        if self.aero is not None:
            self.aero.save_files(case_name=self.case_name, case_route=self.case_route)


def view_wing(case_name, case_route='./', output_folder='./output/'):

    model_settings = {'skin_on': 'off',
                      'discretisation_method': 'michigan',
                      'num_elem': 2}
    pazy = PazyWing(case_name, case_route, in_settings=model_settings)
    pazy.generate_structure()
    pazy.structure.mirror_wing()
    pazy.generate_aero()
    pazy.save_files()

    settings = dict()

    config = configobj.ConfigObj()
    config.filename = './{}.sharpy'.format(case_name)

    settings['SHARPy'] = {
        'flow': ['BeamLoader',
                 'AerogridLoader',
                 'Modal',
                 'BeamPlot',
                 'AerogridPlot'
                 ],
        'case': case_name, 'route': case_route,
        'write_screen': 'on', 'write_log': 'on',
        'log_folder': output_folder + '/' + case_name + '/',
        'log_file': case_name + '.log'}

    settings['BeamLoader'] = {
        'unsteady': 'off'}

    settings['AerogridLoader'] = {
        'unsteady': 'off',
        'aligned_grid': 'on',
        'mstar': 4 * 4,
        'freestream_dir': [1, 0, 0],
        'wake_shape_generator': 'StraightWake',
        'wake_shape_generator_input': {'u_inf': 10,
                                       'u_inf_direction': [1, 0, 0],
                                       'dt': 0.01}}


    settings['Modal'] = {'folder': output_folder,
                         'NumLambda': 20,
                         'rigid_body_modes': 'off',
                         'print_matrices': 'on',
                         'keep_linear_matrices': 'on',
                         'write_dat': 'on',
                         'continuous_eigenvalues': 'off',
                         'dt': 0,
                         'plot_eigenvalues': False,
                         # 'max_rotation_deg': 15.,
                         # 'max_displacement': 0.15,
                         'write_modes_vtk': 'on',
                         'use_undamped_modes': 'on'}

    settings['BeamPlot'] = {'folder': output_folder}

    settings['AerogridPlot'] = {'folder': output_folder,
                              'include_rbm': 'off',
                              'include_applied_forces': 'on',
                              'minus_m_star': 0}

    for k, v in settings.items():
        config[k] = v

    config.write()


if __name__ == '__main__':
    import sharpy.sharpy_main
    # for n_elem in [16, 32, 64, 92, 128]:
    #     pazy = PazyWing()
    #     case_name = 'modal_inertia_even_n{}'.format(n_elem)
    #     output_folder = './output/'
    #     pazy.generate_structure(n_elem, case_name=case_name)
    #     pazy.structure.create_modal_simulation(case_name=case_name, output_folder=output_folder)
    #     sharpy.sharpy_main.main(['', case_name + '.sharpy'])

    # pazy = PazyWing()
    case_name = 'full_wing'
    output_folder = './output/'
    # pazy.generate_structure(n_elem, case_name=case_name)
    # pazy.structure.create_modal_simulation(case_name=case_name, output_folder=output_folder)
    # sharpy.sharpy_main.main(['', case_name + '.sharpy'])

    # model_settings = {'skin_on': 'off',
    #                   'discretisation_method': 'michigan',
    #                   'num_elem': 2}
    #
    # pazy = PazyWing(case_name, case_route='./', in_settings=model_settings)
    # pazy.generate_structure()
    # pazy.structure.mirror_wing()
    # # pazy.structure.add_lumped_mass((10, 0))
    # pazy.save_files()
    # pazy.structure.create_modal_simulation(case_name=case_name, output_folder=output_folder)

    view_wing(case_name, './')
    sharpy.sharpy_main.main(['', case_name + '.sharpy'])

