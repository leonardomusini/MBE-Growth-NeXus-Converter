import os
from pathlib import Path
import h5py
import numpy as np

from parser import *
from utils import *

sample_id_list = file_searcher("./growths2/")

for sample_id in sample_id_list:

    try:
        wri_path = f"./growths2/hm{sample_id}.wri"

        ep_path = None
        for ext in ["ep4", "ep3", "ep2"]:
            ep_path_candidate = f"./growths2/hm{sample_id}.{ext}"
            if os.path.exists(ep_path_candidate):
                ep_path = ep_path_candidate
                break

        refl_path = f"./growths/logs/hm{sample_id}refl.txt"
        pyro_path = f"./growths/logs/hm{sample_id}temp.txt"

        # Required files
        if not os.path.exists(wri_path):
            print("Error: Missing wri file, skipping this iteration.")
            continue  

        if ep_path is None:
            print("Error: Missing ep file, skipping this iteration.")
            continue

        date, starting_time, duration, sample_th, As_pp, Tdeox, rot = parse_wri(wri_path)

        ending_time = time_calculator(starting_time, duration)
        starting_time = time_converter(date, starting_time)
        ending_time = time_converter(date, ending_time)
        h, m, s = map(int, duration.split(":"))
        duration = h + (m / 60) + (s / 3600)

        if isinstance(As_pp, str):
            As_pp = arsenic_ranges(As_pp)

        substrate = parse_substrate(wri_path)

        _, rpm, g_temp, Si_curr, C_curr = parse_wri_layer(wri_path, Tdeox)

        layer, loop, r_step, material, dop_element, thickness, g_time, g_rate, x, shutters, pg_rates = parse_layer(ep_path)

        timestamp_r = wavelengths = dayfraction_r = calibs = None
        timestamp_p = dayfraction_p = t_ratio = t_emiss = None

        # Optional file
        if os.path.exists(refl_path):
            timestamp_r, wavelengths, dayfraction_r, calibs = parse_reflectometer(refl_path)
        #else:
         #  print(f"Warning: Missing {refl_path}, skipping sensor parsing.")
        # Optional file
        if os.path.exists(pyro_path):
            timestamp_p, dayfraction_p, t_ratio, t_emiss = parse_pyrometer(pyro_path)
        #else:
         #  print(f"Warning: Missing {pyro_path}, skipping sensor parsing.")


        # Path of the file
        output_directory = Path("./nexus2")  
        output_directory.mkdir(parents=True, exist_ok=True)

        # File name
        filename = f"hm{sample_id}.nxs"
        file_path = output_directory / filename

        # Creation of the file
        with h5py.File(file_path, "w") as f:

            # NeXus Structure

            # NXentry: /entry
            f.create_group("entry")
            f['/entry'].attrs["NX_class"] = "NXentry"
            f['/entry'].attrs["default"] = "entry"
            # /entry/definition
            f['/entry'].create_dataset('definition',data='NXepitaxy') # application definition
            # /entry/title
            f['/entry'].create_dataset('title', data=f'Growth hm{sample_id} {date}')
            # /entry/experiment_description
            f['/entry'].create_dataset('experiment_description', data='High Mobility Molecular Beam Epitaxy Growth')
            f['/entry/experiment_description'].attrs["description"] = "Growing technique involved"
            # /entry/start_time
            f['/entry'].create_dataset('start_time', data=starting_time)
            # /entry/end_time
            f['/entry'].create_dataset('end_time', data=ending_time)
            # /entry/duration
            f['/entry'].create_dataset('duration', data=duration)
            f['/entry/duration'].attrs["description"] = 'Total time of growth'
            f['/entry/duration'].attrs["units"] = 'hour'


            date = datetime.strptime(date, "%m-%d-%Y")
            if date < datetime(2023,11,2):

                # NXuser: /entry/user
                f['/entry'].create_group("user")
                f['/entry/user'].attrs["NX_class"] = "NXuser"
                # /entry/user/name
                f['/entry/user'].create_dataset('name',data='Giorgio Biasiol')
                # /entry/user/affiliation
                f['/entry/user'].create_dataset('affiliation',data='CNR-IOM')
                # /entry/user/role
                f['/entry/user'].create_dataset('role',data='Technology Director')
                # /entry/user/email
                f['/entry/user'].create_dataset('email',data='giorgio.biasiol@cnr.it')
                # /entry/user/ORCID
                f['/entry/user'].create_dataset('ORCID',data='https://orcid.org/0000-0001-7974-5459')

            else:

                # NXuser: /entry/user_1
                f['/entry'].create_group("user_1")
                f['/entry/user_1'].attrs["NX_class"] = "NXuser"
                # /entry/user/name
                f['/entry/user_1'].create_dataset('name',data='Giorgio Biasiol')
                # /entry/user/affiliation
                f['/entry/user_1'].create_dataset('affiliation',data='CNR-IOM')
                # /entry/user/role
                f['/entry/user_1'].create_dataset('role',data='Technology Director')
                # /entry/user/email
                f['/entry/user_1'].create_dataset('email',data='giorgio.biasiol@cnr.it')
                # /entry/user/ORCID
                f['/entry/user_1'].create_dataset('ORCID',data='https://orcid.org/0000-0001-7974-5459')

                # NXuser: /entry/user_2
                f['/entry'].create_group("user_2")
                f['/entry/user_2'].attrs["NX_class"] = "NXuser"
                # /entry/user/name
                f['/entry/user_2'].create_dataset('name',data='Davide Curcio')
                # /entry/user/affiliation
                f['/entry/user_2'].create_dataset('affiliation',data='CNR-IOM')
                # /entry/user/role
                f['/entry/user_2'].create_dataset('role',data='Research technologist')
                # /entry/user/email
                f['/entry/user_2'].create_dataset('email',data='curcio@iom.cnr.it')
                # /entry/user/ORCID
                f['/entry/user_2'].create_dataset('ORCID',data='https://orcid.org/0000-0003-2488-3840')


            # NXsample: /entry/sample
            f['/entry'].create_group("sample")
            f['/entry/sample'].attrs["NX_class"] = "NXsample"
            # /entry/sample/name
            f['/entry/sample'].create_dataset('name', data=f'hm{sample_id}')
            # /entry/sample/thickness
            f['/entry/sample'].create_dataset('thickness', data=sample_th)
            f['/entry/sample/thickness'].attrs['units'] = 'μm'

            # NXsample: /entry/sample/substrate
            f['/entry/sample'].create_group("substrate")
            f['/entry/sample/substrate'].attrs["NX_class"] = "NXsample_component"
            # /entry/sample/substrate/name
            f['/entry/sample/substrate'].create_dataset('name', data=substrate['name'])
            f['/entry/sample/substrate/name'].attrs["description"] = 'Wafer model'
            # /entry/sample/substrate/chemical_formula
            f['/entry/sample/substrate'].create_dataset('chemical_formula', data='GaAs')
            # /entry/sample/substrate/thickness
            if np.isnan(substrate['thickness']):
                f['/entry/sample/substrate'].create_dataset('thickness', data=500)
            else:
                f['/entry/sample/substrate'].create_dataset('thickness', data=substrate['thickness'])
            f['/entry/sample/substrate/thickness'].attrs["description"] = 'Thickness of the substrate'
            f['/entry/sample/substrate/thickness'].attrs['units'] = 'μm'
            # /entry/sample/substrate/area
            f['/entry/sample/substrate'].create_dataset('area', data=area_converter(substrate['area']))
            f['/entry/sample/substrate/area'].attrs['description'] = 'Area of the substrate'
            f['/entry/sample/substrate/area'].attrs['units'] = 'mm^2'
            # /entry/sample/substrate/diameter
            f['/entry/sample/substrate'].create_dataset('diameter', data=substrate['diameter'])
            f['/entry/sample/substrate/diameter'].attrs['description'] = 'Wafer diameter'
            f['/entry/sample/substrate/diameter'].attrs['units'] = 'in'
            # /entry/sample/substrate/crystal_orientation
            f['/entry/sample/substrate'].create_dataset('crystal_orientation', data=substrate['orientation'])
            f['/entry/sample/substrate/crystal_orientation'].attrs['description'] = 'Crystallographic direction of the material'
            # /entry/sample/substrate/flat_convention
            f['/entry/sample/substrate'].create_dataset('flat_convention', data=substrate['convention'])
            f['/entry/sample/substrate/flat_convention'].attrs['description'] = 'Flat convention of the wafer, value accepted are EJ or US' 
            # /entry/sample/substrate/doping
            f['/entry/sample/substrate'].create_dataset('doping', data=substrate['doping'])
            f['/entry/sample/substrate/doping'].attrs['description'] = 'Doping type and level of the substrate (e.g., p+, n, SI)'
            # /entry/sample/substrate/holder
            f['/entry/sample/substrate'].create_dataset('holder', data=substrate['holder'])
            f['/entry/sample/substrate/holder'].attrs['description'] = 'Type of substrate holder used in the process' 
            # /entry/sample/substrate/crystalline_structure
            f['/entry/sample/substrate'].create_dataset('crystalline_structure', data='single crystal')
            f['/entry/sample/substrate/crystalline_structure'].attrs['description'] = 'Crystalline structure of the material'

            # Initialize the current layer counter
            current_layer = 0

            # Iterate through each layer
            for i in range(len(layer)):
                layer_name = f"layer{current_layer + 1:02d}"
                 # NXsample_component: /entry/sample/layer
                f['/entry/sample/'].create_group(f"{layer_name}")
                f[f'/entry/sample/{layer_name}'].attrs["NX_class"] = "NXsample_component"
                # /entry/sample/layer/name
                f[f'/entry/sample/{layer_name}'].create_dataset('name', data=layer_name)
                f[f'/entry/sample/{layer_name}/name'].attrs["description"] = 'Name of the layer'
                # /entry/sample/layer/chemical_formula
                f[f'/entry/sample/{layer_name}'].create_dataset('chemical_formula', data=material[i])
                f[f'/entry/sample/{layer_name}/chemical_formula'].attrs["description"] = 'Chemical formula of the layer'
                # /entry/sample/layer/doping
                if material[i] == "Si":
                    f[f'/entry/sample/{layer_name}'].create_dataset('doping', data=g_time[i]*10e11)
                    f[f'/entry/sample/{layer_name}/doping'].attrs["description"] = 'Doping value of the layer'
                    f[f'/entry/sample/{layer_name}/doping'].attrs["units"] = 'cm^-2'
                else:
                    f[f'/entry/sample/{layer_name}'].create_dataset('doping', data=doping_calculator(dop_element[i],Si_curr[i],C_curr[i],g_rate[i]))
                    f[f'/entry/sample/{layer_name}/doping'].attrs["description"] = 'Doping value of the layer'
                    f[f'/entry/sample/{layer_name}/doping'].attrs["units"] = 'cm^-3'
                # /entry/sample/layer/thickness
                f[f'/entry/sample/{layer_name}'].create_dataset('thickness', data=thickness[i])
                f[f'/entry/sample/{layer_name}/thickness'].attrs["units"] = 'Å'
                # /entry/sample/layer/growth_temperature
                f[f'/entry/sample/{layer_name}'].create_dataset('growth_temperature', data=g_temp[i])
                f[f'/entry/sample/{layer_name}/growth_temperature'].attrs["description"] = 'Growing temperature of the layer'
                f[f'/entry/sample/{layer_name}/growth_temperature'].attrs["units"] = '°C'
                # /entry/sample/layer/growth_time
                f[f'/entry/sample/{layer_name}'].create_dataset('growth_time', data=g_time[i])
                f[f'/entry/sample/{layer_name}/growth_time'].attrs["description"] = 'Growing time of the layer'
                f[f'/entry/sample/{layer_name}/growth_time'].attrs["units"] = 's'
                # /entry/sample/layer/growth_rate
                f[f'/entry/sample/{layer_name}'].create_dataset('growth_rate', data=g_rate[i])
                f[f'/entry/sample/{layer_name}/growth_rate'].attrs["description"] = 'Rate at which the layer is grown'
                f[f'/entry/sample/{layer_name}/growth_rate'].attrs["units"] = 'Å/s'
                # /entry/sample/layer/alloy_fraction
                f[f'/entry/sample/{layer_name}'].create_dataset('alloy_fraction', data=x[i])
                f[f'/entry/sample/{layer_name}/alloy_fraction'].attrs["description"] = 'Fraction of the first element in a ternary alloy'
                # /entry/sample/layer/rotational_frequency
                if sum(rpm) == 0 and not np.isnan(rot):
                    f[f'/entry/sample/{layer_name}'].create_dataset('rotational_frequency', data=rot)
                else:
                    f[f'/entry/sample/{layer_name}'].create_dataset('rotational_frequency', data=rpm[i])
                f[f'/entry/sample/{layer_name}/rotational_frequency'].attrs["description"] = 'Rotational frequency of the sample for the current layer'
                f[f'/entry/sample/{layer_name}/rotational_frequency'].attrs["units"] = 'rpm'

                cell_material = ["Ga", "Ga", "Al", "In"]
                cell_names = ["Cell Ga1", "Cell Ga2", "Cell Al", "Cell In"]
                cell_number = 0
                for k, open in enumerate(shutters[i]):
                    if open:
                        cell_number += 1

                        # NXmaterial_source: /entry/sample/layer/material_source
                        f[f'/entry/sample/{layer_name}'].create_group(f"cell_{cell_number}")
                        f[f'/entry/sample/{layer_name}/cell_{cell_number}'].attrs["NX_class"] = "NXmaterial_source"
                        # /entry/sample/layer/cell/name
                        f[f'/entry/sample/{layer_name}/cell_{cell_number}'].create_dataset('name', data=cell_names[k])
                        f[f'/entry/sample/{layer_name}/cell_{cell_number}/name'].attrs["description"] = 'Name of the material source device'
                        # /entry/sample/layer/cell/model
                        f[f'/entry/sample/{layer_name}/cell_{cell_number}'].create_dataset('model', data='')
                        f[f'/entry/sample/{layer_name}/cell_{cell_number}/model'].attrs["description"] = 'Model of the material source device'
                        # /entry/sample/layer/cell/type
                        f[f'/entry/sample/{layer_name}/cell_{cell_number}'].create_dataset('type', data='effusion_cell')
                        f[f'/entry/sample/{layer_name}/cell_{cell_number}/type'].attrs["description"] = 'Type of material source device'
                        # /entry/sample/layer/cell/shutter_status
                        f[f'/entry/sample/{layer_name}/cell_{cell_number}'].create_dataset('shutter_status', data='open')
                        f[f'/entry/sample/{layer_name}/cell_{cell_number}/shutter_status'].attrs["description"] = 'Status of the shutter during this layer growth step'
                        # /entry/sample/layer/cell/material
                        f[f'/entry/sample/{layer_name}/cell_{cell_number}'].create_dataset('material', data=cell_material[k])
                        f[f'/entry/sample/{layer_name}/cell_{cell_number}/material'].attrs["description"] = 'Material of the cell'
                        # /entry/sample/layer/cell/partial_growth_rate
                        f[f'/entry/sample/{layer_name}/cell_{cell_number}'].create_dataset('partial_growth_rate', data=pg_rates[i][k])
                        f[f'/entry/sample/{layer_name}/cell_{cell_number}'].attrs["description"] = 'Partial growth rate of the cell'
                        f[f'/entry/sample/{layer_name}/cell_{cell_number}'].attrs["units"] = 'Å/s'

                # NXmaterial_source: /entry/sample/layer/material_source
                f[f'/entry/sample/{layer_name}'].create_group(f"cell_{cell_number+1}")
                f[f'/entry/sample/{layer_name}/cell_{cell_number+1}'].attrs["NX_class"] = "NXmaterial_source"
                # /entry/sample/layer/cell/name
                f[f'/entry/sample/{layer_name}/cell_{cell_number+1}'].create_dataset('name', data='Cell As')
                f[f'/entry/sample/{layer_name}/cell_{cell_number+1}/name'].attrs["description"] = 'Name of the material source device'
                # /entry/sample/layer/cell/model
                f[f'/entry/sample/{layer_name}/cell_{cell_number+1}'].create_dataset('model', data='')
                f[f'/entry/sample/{layer_name}/cell_{cell_number+1}/model'].attrs["description"] = 'Model of the material source device'
                # /entry/sample/layer/cell/type
                f[f'/entry/sample/{layer_name}/cell_{cell_number+1}'].create_dataset('type', data='effusion_cell')
                f[f'/entry/sample/{layer_name}/cell_{cell_number+1}/type'].attrs["description"] = 'Type of material source device'
                # /entry/sample/layer/cell/shutter_status
                f[f'/entry/sample/{layer_name}/cell_{cell_number+1}'].create_dataset('shutter_status', data='open')
                f[f'/entry/sample/{layer_name}/cell_{cell_number+1}/shutter_status'].attrs["description"] = 'Status of the shutter during this layer growth step'
                # /entry/sample/layer/cell/material
                f[f'/entry/sample/{layer_name}/cell_{cell_number+1}'].create_dataset('material', data='As')
                f[f'/entry/sample/{layer_name}/cell_{cell_number+1}/material'].attrs["description"] = 'Material of the cell'
                # /entry/sample/layer/cell/partial_pressure
                f[f'/entry/sample/{layer_name}/cell_{cell_number+1}'].create_dataset('partial_growth_rate', data=As_pp)
                f[f'/entry/sample/{layer_name}/cell_{cell_number+1}'].attrs["description"] = 'Partial pressure of the cell'
                f[f'/entry/sample/{layer_name}/cell_{cell_number+1}'].attrs["units"] = 'Torr'

                current_layer += 1

                # --- Handle Repetition ---

                if loop[i] > 1 and i > 0:  
                    repeat_count = loop[i] - 1  # Until now there is already an occurence of the layers to repeat
                    for repeat in range(repeat_count):
                        # Repeat layers starting from the return step r_step[i]
                        for j in range(r_step[i]-1, i+1):
                            layer_name = f"layer{current_layer + 1:02d}"
                            # NXsample_component: /entry/sample/layer
                            f['/entry/sample/'].create_group(f"{layer_name}")
                            f[f'/entry/sample/{layer_name}'].attrs["NX_class"] = "NXsample_component"
                            # /entry/sample/layer/name
                            f[f'/entry/sample/{layer_name}'].create_dataset('name', data=layer_name)
                            # /entry/sample/layer/chemical_formula
                            f[f'/entry/sample/{layer_name}'].create_dataset('chemical_formula', data=material[j])
                            if material[j] == "Si":
                                f[f'/entry/sample/{layer_name}'].create_dataset('doping', data=g_time[j]*10e11)
                                f[f'/entry/sample/{layer_name}/doping'].attrs["description"] = 'Doping value of the layer'
                                f[f'/entry/sample/{layer_name}/doping'].attrs["units"] = 'cm^-2'
                            else:
                                f[f'/entry/sample/{layer_name}'].create_dataset('doping', data=doping_calculator(dop_element[j],Si_curr[j],C_curr[j],g_rate[j]))
                                f[f'/entry/sample/{layer_name}/doping'].attrs["description"] = 'Doping value of the layer'
                                f[f'/entry/sample/{layer_name}/doping'].attrs["units"] = 'cm^-3'
                            # /entry/sample/layer/thickness
                            f[f'/entry/sample/{layer_name}'].create_dataset('thickness', data=thickness[j])
                            f[f'/entry/sample/{layer_name}/thickness'].attrs["units"] = 'Å'
                            # /entry/sample/layer/growth_temperature
                            f[f'/entry/sample/{layer_name}'].create_dataset('growth_temperature', data=g_temp[j])
                            f[f'/entry/sample/{layer_name}/growth_temperature'].attrs["description"] = 'Growing temperature of the layer'
                            f[f'/entry/sample/{layer_name}/growth_temperature'].attrs["units"] = '°C'
                            # /entry/sample/layer/growth_time
                            f[f'/entry/sample/{layer_name}'].create_dataset('growth_time', data=g_time[j])
                            f[f'/entry/sample/{layer_name}/growth_time'].attrs["description"] = 'Growing time of the layer'
                            f[f'/entry/sample/{layer_name}/growth_time'].attrs["units"] = 's'
                            # /entry/sample/layer/growth_rate
                            f[f'/entry/sample/{layer_name}'].create_dataset('growth_rate', data=g_rate[j])
                            f[f'/entry/sample/{layer_name}/growth_rate'].attrs["description"] = 'Rate at which the layer is grown'
                            f[f'/entry/sample/{layer_name}/growth_rate'].attrs["units"] = 'Å/s'
                            # /entry/sample/layer/alloy_fraction
                            f[f'/entry/sample/{layer_name}'].create_dataset('alloy_fraction', data=x[j])
                            f[f'/entry/sample/{layer_name}/alloy_fraction'].attrs["description"] = 'Fraction of the first element in a ternary alloy'
                            # /entry/sample/layer/rotational_frequency
                            if sum(rpm) == 0 and not np.isnan(rot):
                                f[f'/entry/sample/{layer_name}'].create_dataset('rotational_frequency', data=rot)
                            else:
                                f[f'/entry/sample/{layer_name}'].create_dataset('rotational_frequency', data=rpm[j])
                            f[f'/entry/sample/{layer_name}/rotational_frequency'].attrs["description"] = 'Rotational frequency of the sample for the current layer'
                            f[f'/entry/sample/{layer_name}/rotational_frequency'].attrs["units"] = 'rpm'
                            
                            cell_material = ["Ga", "Ga", "Al", "In"]
                            cell_names = ["Cell Ga1", "Cell Ga2", "Cell Al", "Cell In"]
                            cell_number = 0
                            for k, open in enumerate(shutters[i]):
                                if open:
                                    cell_number += 1

                                    # NXmaterial_source: /entry/sample/layer/material_source
                                    f[f'/entry/sample/{layer_name}'].create_group(f"cell_{cell_number}")
                                    f[f'/entry/sample/{layer_name}/cell_{cell_number}'].attrs["NX_class"] = "NXmaterial_source"
                                    # /entry/sample/layer/cell/name
                                    f[f'/entry/sample/{layer_name}/cell_{cell_number}'].create_dataset('name', data=cell_names[k])
                                    f[f'/entry/sample/{layer_name}/cell_{cell_number}/name'].attrs["description"] = 'Name of the material source device'
                                    # /entry/sample/layer/cell/model
                                    f[f'/entry/sample/{layer_name}/cell_{cell_number}'].create_dataset('model', data='')
                                    f[f'/entry/sample/{layer_name}/cell_{cell_number}/model'].attrs["description"] = 'Model of the material source device'
                                    # /entry/sample/layer/cell/type
                                    f[f'/entry/sample/{layer_name}/cell_{cell_number}'].create_dataset('type', data='effusion_cell')
                                    f[f'/entry/sample/{layer_name}/cell_{cell_number}/type'].attrs["description"] = 'Type of material source device'
                                    # /entry/sample/layer/cell/shutter_status
                                    f[f'/entry/sample/{layer_name}/cell_{cell_number}'].create_dataset('shutter_status', data='open')
                                    f[f'/entry/sample/{layer_name}/cell_{cell_number}/shutter_status'].attrs["description"] = 'Status of the shutter during this layer growth step'
                                    # /entry/sample/layer/cell/material
                                    f[f'/entry/sample/{layer_name}/cell_{cell_number}'].create_dataset('material', data=cell_material[k])
                                    f[f'/entry/sample/{layer_name}/cell_{cell_number}/material'].attrs["description"] = 'Material of the cell'
                                    # /entry/sample/layer/cell/partial_growth_rate
                                    f[f'/entry/sample/{layer_name}/cell_{cell_number}'].create_dataset('partial_growth_rate', data=pg_rates[i][k])
                                    f[f'/entry/sample/{layer_name}/cell_{cell_number}'].attrs["description"] = 'Partial growth rate of the cell'
                                    f[f'/entry/sample/{layer_name}/cell_{cell_number}'].attrs["units"] = 'Å/s'

                            # NXmaterial_source: /entry/sample/layer/material_source
                            f[f'/entry/sample/{layer_name}'].create_group(f"cell_{cell_number+1}")
                            f[f'/entry/sample/{layer_name}/cell_{cell_number+1}'].attrs["NX_class"] = "NXmaterial_source"
                            # /entry/sample/layer/cell/name
                            f[f'/entry/sample/{layer_name}/cell_{cell_number+1}'].create_dataset('name', data='Cell As')
                            f[f'/entry/sample/{layer_name}/cell_{cell_number+1}/name'].attrs["description"] = 'Name of the material source device'
                            # /entry/sample/layer/cell/model
                            f[f'/entry/sample/{layer_name}/cell_{cell_number+1}'].create_dataset('model', data='')
                            f[f'/entry/sample/{layer_name}/cell_{cell_number+1}/model'].attrs["description"] = 'Model of the material source device'
                            # /entry/sample/layer/cell/type
                            f[f'/entry/sample/{layer_name}/cell_{cell_number+1}'].create_dataset('type', data='effusion_cell')
                            f[f'/entry/sample/{layer_name}/cell_{cell_number+1}/type'].attrs["description"] = 'Type of material source device'
                            # /entry/sample/layer/cell/shutter_status
                            f[f'/entry/sample/{layer_name}/cell_{cell_number+1}'].create_dataset('shutter_status', data='open')
                            f[f'/entry/sample/{layer_name}/cell_{cell_number+1}/shutter_status'].attrs["description"] = 'Status of the shutter during this layer growth step'
                            # /entry/sample/layer/cell/material
                            f[f'/entry/sample/{layer_name}/cell_{cell_number+1}'].create_dataset('material', data='As')
                            f[f'/entry/sample/{layer_name}/cell_{cell_number+1}/material'].attrs["description"] = 'Material of the cell'
                            # /entry/sample/layer/cell/partial_pressure
                            f[f'/entry/sample/{layer_name}/cell_{cell_number+1}'].create_dataset('partial_growth_rate', data=As_pp)
                            f[f'/entry/sample/{layer_name}/cell_{cell_number+1}'].attrs["description"] = 'Partial pressure of the cell'
                            f[f'/entry/sample/{layer_name}/cell_{cell_number+1}'].attrs["units"] = 'Torr'

                            current_layer += 1

            

            # NXinstrument: /entry/instrument
            f['/entry'].create_group("instrument")
            f['/entry/instrument'].attrs["NX_class"] = "NXinstrument"

            # /entry/instrument/chamber
            f['/entry/instrument'].create_group("chamber")
            f['/entry/instrument/chamber'].attrs["NX_class"] = "NXenvironment"
            # /entry/instrument/chamber/name
            f['/entry/instrument/chamber'].create_dataset('name',data='Built in-house')
            f['/entry/instrument/chamber/name'].attrs['description'] = 'Model of the growing chamber'
            # /entry/instrument/chamber/type
            f['/entry/instrument/chamber'].create_dataset('type',data='Ultra High Vacuum (UHV) Chamber')
            f['/entry/instrument/chamber/type'].attrs['description'] = 'Type of growing chamber'
            # /entry/instrument/chamber/description
            f['/entry/instrument/chamber'].create_dataset('description',data='')
            f['/entry/instrument/chamber/description'].attrs['description'] = 'Information of the chamber (ex. cooling system)'

            # /entry/instrument/chamber/sensor_1
            f['/entry/instrument/chamber'].create_group(f"sensor_1")
            f[f'/entry/instrument/chamber/sensor_1'].attrs["NX_class"] = "NXsensor"
            # /entry/instrument/chamber/sensor_1/name
            f[f'/entry/instrument/chamber/sensor_1'].create_dataset('name', data='Reflectometer 950')
            f[f'/entry/instrument/chamber/sensor_1/name'].attrs['description'] = 'Name of sensor'
            # /entry/instrument/chamber/sensor_1/model
            f[f'/entry/instrument/chamber/sensor_1'].create_dataset('model', data='')
            f[f'/entry/instrument/chamber/sensor_1/model'].attrs['description'] = 'Model of sensor'
            # /entry/instrument/chamber/sensor_1/measurement
            f[f'/entry/instrument/chamber/sensor_1'].create_dataset('measurement', data='reflectivity')
            f[f'/entry/instrument/chamber/sensor_1/measurement'].attrs['description'] = 'Physical quantity being measured'
            if os.path.exists(refl_path):
                # /entry/instrument/chamber/sensor_1/reflectivity
                f[f'/entry/instrument/chamber/sensor_1'].create_group("reflectivity")
                f[f'/entry/instrument/chamber/sensor_1/reflectivity'].attrs["NX_class"] = "NXlog"
                # /entry/instrument/chamber/sensor_1/reflectivity/time
                f[f'/entry/instrument/chamber/sensor_1/reflectivity'].create_dataset('time', data=dayfraction_r)
                f[f'/entry/instrument/chamber/sensor_1/reflectivity/time'].attrs['start'] = timestamp_r
                f[f'/entry/instrument/chamber/sensor_1/reflectivity/time'].attrs['unit'] = 'day fraction'
                # /entry/instrument/chamber/sensor_1/reflectivity/value
                f[f'/entry/instrument/chamber/sensor_1/reflectivity'].create_dataset('value', data=calibs[0])
                f[f'/entry/instrument/chamber/sensor_1/reflectivity/value'].attrs['description'] = 'Reflectivity of the sample'

            # /entry/instrument/chamber/sensor_2
            f['/entry/instrument/chamber'].create_group(f"sensor_2")
            f[f'/entry/instrument/chamber/sensor_2'].attrs["NX_class"] = "NXsensor"
            # /entry/instrument/chamber/sensor_2/name
            f[f'/entry/instrument/chamber/sensor_2'].create_dataset('name', data='Reflectometer 470')
            f[f'/entry/instrument/chamber/sensor_2/name'].attrs['description'] = 'Name of sensor'
            # /entry/instrument/chamber/sensor_2/model
            f[f'/entry/instrument/chamber/sensor_2'].create_dataset('model', data='')
            f[f'/entry/instrument/chamber/sensor_2/model'].attrs['description'] = 'Model of sensor'
            # /entry/instrument/chamber/sensor_2/measurement
            f[f'/entry/instrument/chamber/sensor_2'].create_dataset('measurement', data='reflectivity')
            f[f'/entry/instrument/chamber/sensor_2/measurement'].attrs['description'] = 'Physical quantity being measured'
            if os.path.exists(refl_path):
                # /entry/instrument/chamber/sensor_2/reflectivity
                f[f'/entry/instrument/chamber/sensor_2'].create_group("reflectivity")
                f[f'/entry/instrument/chamber/sensor_2/reflectivity'].attrs["NX_class"] = "NXlog"
                # /entry/instrument/chamber/sensor_2/reflectivity/time
                f[f'/entry/instrument/chamber/sensor_2/reflectivity'].create_dataset('time', data=dayfraction_r)
                f[f'/entry/instrument/chamber/sensor_2/reflectivity/time'].attrs['start'] = timestamp_r
                f[f'/entry/instrument/chamber/sensor_2/reflectivity/time'].attrs['unit'] = 'day fraction'
                # /entry/instrument/chamber/sensor_2/reflectivity/value
                f[f'/entry/instrument/chamber/sensor_2/reflectivity'].create_dataset('value', data=calibs[1])
                f[f'/entry/instrument/chamber/sensor_2/reflectivity/value'].attrs['description'] = 'Reflectivity of the sample'

            # /entry/instrument/chamber/sensor_3
            f['/entry/instrument/chamber'].create_group("sensor_3")
            f['/entry/instrument/chamber/sensor_3'].attrs["NX_class"] = "NXsensor"
            # /entry/instrument/chamber/sensor_3/name
            f['/entry/instrument/chamber/sensor_3'].create_dataset('name', data='One-channel pyrometer')
            f['/entry/instrument/chamber/sensor_3/name'].attrs['description'] = 'Name of the sensor'
            # /entry/instrument/chamber/sensor_3/model
            f['/entry/instrument/chamber/sensor_3'].create_dataset('model', data='')
            f['/entry/instrument/chamber/sensor_3/model'].attrs['description'] = 'Model of the sensor'
            # /entry/instrument/chamber/sensor_3/measurement
            f['/entry/instrument/chamber/sensor_3'].create_dataset('measurement', data='rate_temperature')
            f['/entry/instrument/chamber/sensor_3/measurement'].attrs['description'] = 'Physical quantity being measured'
            if os.path.exists(pyro_path):
                # /entry/instrument/chamber/sensor_3/temperature
                f['/entry/instrument/chamber/sensor_3'].create_group("temperature")
                f['/entry/instrument/chamber/sensor_3/temperature'].attrs["NX_class"] = "NXlog"
                # /entry/instrument/chamber/sensor_3/temperature/time
                f['/entry/instrument/chamber/sensor_3/temperature'].create_dataset('time', data=dayfraction_p)
                f['/entry/instrument/chamber/sensor_3/temperature/time'].attrs['start'] = timestamp_p
                f['/entry/instrument/chamber/sensor_3/temperature/time'].attrs['unit'] = 'day fraction'
                # /entry/instrument/chamber/sensor_3/temperature/value
                f['/entry/instrument/chamber/sensor_3/temperature'].create_dataset('value', data=t_ratio)
                f['/entry/instrument/chamber/sensor_3/temperature/value'].attrs['description'] = 'Rate temperature'

            # /entry/instrument/chamber/sensor_4
            f['/entry/instrument/chamber'].create_group("sensor_4")
            f['/entry/instrument/chamber/sensor_4'].attrs["NX_class"] = "NXsensor"
            # /entry/instrument/chamber/sensor_4/name
            f['/entry/instrument/chamber/sensor_4'].create_dataset('name', data='Two-channel pyrometer')
            f['/entry/instrument/chamber/sensor_4/name'].attrs['description'] = 'Name of the sensor'
            # /entry/instrument/chamber/sensor_4/model
            f['/entry/instrument/chamber/sensor_4'].create_dataset('model', data='')
            f['/entry/instrument/chamber/sensor_4/model'].attrs['description'] = 'Model of the sensor'
            # /entry/instrument/chamber/sensor_4/measurement
            f['/entry/instrument/chamber/sensor_4'].create_dataset('measurement', data='emissivity_temperature')
            f['/entry/instrument/chamber/sensor_4/measurement'].attrs['description'] = 'Physical quantity being measured'
            if os.path.exists(pyro_path):
                # /entry/instrument/chamber/sensor_4/temperature
                f['/entry/instrument/chamber/sensor_4'].create_group("temperature")
                f['/entry/instrument/chamber/sensor_4/temperature'].attrs["NX_class"] = "NXlog"
                # /entry/instrument/chamber/sensor_4/temperature/time
                f['/entry/instrument/chamber/sensor_4/temperature'].create_dataset('time', data=dayfraction_p)
                f['/entry/instrument/chamber/sensor_4/temperature/time'].attrs['start'] = timestamp_p
                f['/entry/instrument/chamber/sensor_2/temperature/time'].attrs['unit'] = 'day fraction'
                # /entry/instrument/chamber/sensor_4/temperature/value
                f['/entry/instrument/chamber/sensor_4/temperature'].create_dataset('value', data=t_emiss)
                f['/entry/instrument/chamber/sensor_4/temperature/value'].attrs['description'] = 'Emissivity temperature'
                f['/entry/instrument/chamber/sensor_4/temperature/value'].attrs['unit'] = '°C'

            # /entry/instrument/chamber/sensor_5
            f['/entry/instrument/chamber'].create_group("sensor_5")
            f['/entry/instrument/chamber/sensor_5'].attrs['NX_class'] = "NXsensor" 
            # /entry/instrument/chamber/sensor_5/name
            f['/entry/instrument/chamber/sensor_5'].create_dataset('name',data='Ion gauge')
            f['/entry/instrument/chamber/sensor_5/name'].attrs['description'] = 'Name of the ion gauge'
            # /entry/instrument/chamber/sensor_5/model
            f['/entry/instrument/chamber/sensor_5'].create_dataset('model',data='')
            f['/entry/instrument/chamber/sensor_5/model'].attrs['description'] = 'Model of the ion gauge' 
            # /entry/instrument/chamber/sensor_5/measurement
            f['/entry/instrument/chamber/sensor_5'].create_dataset('measurement',data='pressure')
            f['/entry/instrument/chamber/sensor_5/measurement'].attrs['description'] = 'Physical quantity being measured'
            # /entry/instrument/chamber/sensor_5/value
            f['/entry/instrument/chamber/sensor_5'].create_dataset('value',data=1e-12)
            f['/entry/instrument/chamber/sensor_5/value'].attrs['description'] = 'Nominal value of the signal'
            f['/entry/instrument/chamber/sensor_5/value'].attrs['unit'] = 'mbar'

            # NXcooling_device: /entry/instrument/cooling_device
            f['/entry/instrument'].create_group("cooling_device")
            f['/entry/instrument/cooling_device'].attrs["NX_class"] = "NXcooling_device"
            # /entry/instrument/cooling_device/name
            f['/entry/instrument/cooling_device'].create_dataset('name',data='Cooling shroud')
            f['/entry/instrument/cooling_device/name'].attrs['description'] = 'Name of the cooling system'
            # /entry/instrument/cooling_device/model
            f['/entry/instrument/cooling_device'].create_dataset('model',data='')
            f['/entry/instrument/cooling_device/model'].attrs['description'] = 'Model of the cooling system'
            # /entry/instrument/cooling_device/cooling_mode
            f['/entry/instrument/cooling_device'].create_dataset('cooling_mode',data='liquid_nitrogen')
            f['/entry/instrument/cooling_device/cooling_mode'].attrs['description'] = 'Mode used to cool the chamber'
            # /entry/instrument/cooling_device/temperature
            f['/entry/instrument/cooling_device'].create_dataset('temperature',data=77)
            f['/entry/instrument/cooling_device/temperature'].attrs['description'] = 'Nominal temperature of the system'
            f['/entry/instrument/cooling_device/temperature'].attrs['units'] = 'K'

        
        if os.path.exists(file_path):
            print(f"Success: File {filename} created successfully at {file_path}.")
        else:
            print(f"Error: Failed to create {filename}.")

    except Exception as e:
        print(f"Error processing hm{sample_id}: {e}")