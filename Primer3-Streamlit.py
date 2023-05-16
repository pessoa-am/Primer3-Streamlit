# Copyright 2023 Alberto Pessoa <pessoa.am@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# # You should have received a copy of the
# GNU Affero General Public License along with this program.
# If not, see <https://www.gnu.org/licenses/>.


import streamlit as st
from primer3.bindings import design_primers
from primer3_st_args import example_values, task_help, primer_task, table_th, table_salt, st_args, st_default_values, st_static_values
import misprime_libs

colors = {'LEFT': '#8383FC', 'INTERNAL': '#DE83FC', 'RIGHT': '#FCFC83', 'EXCLUDED_REGION': '#FB6A6A', 'TARGET': '#0DF20D', 'INCLUDED_REGION': '#83FCFC'}


def hierarchize(primer3_results):
    results = {}
    p = 0
    results['PRIMERS'] = []
    results['EXPLAIN'] = {}
    for k, v in primer3_results.items():
        if k.startswith('PRIMER_'):
            tmp = k.split('_', 3)
            if len(tmp) == 2:
                results[tmp[1]] = v
                continue
            if tmp[1] == 'RIGHTINTERNAL':
                tmp[1] = 'INTERNAL'
            if tmp[2] == 'EXPLAIN':
                results['EXPLAIN'][tmp[0] + '_' + tmp[1]] = v
            elif tmp[2].isdigit():
                p = int(tmp[2])
                while p + 1 > len(results['PRIMERS']):
                    results['PRIMERS'].append({})
                if tmp[1] not in results['PRIMERS'][p]:
                    results['PRIMERS'][p][tmp[1]] = {}
                if len(tmp) > 3:
                    results['PRIMERS'][p][tmp[1]][tmp[3]] = v
                else:
                    results['PRIMERS'][p][tmp[1]]['POSITION'] = v[0]
                    results['PRIMERS'][p][tmp[1]]['LENGTH'] = v[1]
            elif '{' in tmp[2]:
                results['PRIMERS'][p][tmp[1]][tmp[3]] = v
            elif tmp[2] == 'NUM':
                pass
        else:
            results[k] = v
    for p in results['PRIMERS']:
        if 'LEFT' in p:
            if len(p['LEFT'].keys()) <= 1:
                del p['LEFT']
        if 'RIGHT' in p:
            if len(p['RIGHT'].keys()) <= 1:
                del p['RIGHT']
    return results


def sequence_block(seq, LEFT_POSITION=None, LEFT_LENGTH=None, INTERNAL_POSITION=None, INTERNAL_LENGTH=None, RIGHT_POSITION=None, RIGHT_LENGTH=None, EXCLUDED_REGION=None, TARGET=None, INCLUDED_REGION=None):
    def non_intersecting_parts(ranges, ref):
        non_intersecting = []
        label, start, end = ref
        for _, r_start, r_end in ranges:
            if r_end >= start and r_start <= end:
                if r_start > start:
                    non_intersecting.append((label, start, r_start - 1))
                start = r_end + 1
        if end >= start:
            non_intersecting.append((label, start, end))
        return non_intersecting

    def remove_overlap(ranges):
        ranges.sort(key=lambda x: x[1])
        merged = [ranges[0]]
        for current in ranges[1:]:
            previous = merged[-1]
            if current[1] <= previous[2]:
                merged[-1] = (previous[0], previous[1], max(previous[2], current[2]))
            else:
                merged.append(current)
        return merged

    ranges = []
    io_range = None
    p_ranges = []
    e_ranges = []
    t_ranges = []
    range_gaps = []
    block_coords = []
    ir_start = len(seq) + 1
    ir_end = 0
    if LEFT_POSITION is not None and LEFT_LENGTH is not None:
        p_ranges.append(('LEFT', LEFT_POSITION - 1, LEFT_POSITION - 1 + LEFT_LENGTH - 1))
    if INTERNAL_POSITION is not None and INTERNAL_LENGTH is not None:
        io_range = ('INTERNAL', INTERNAL_POSITION - 1, INTERNAL_POSITION - 1 + INTERNAL_LENGTH - 1)
    if RIGHT_POSITION is not None and RIGHT_LENGTH is not None:
        p_ranges.append(('RIGHT', RIGHT_POSITION - RIGHT_LENGTH, RIGHT_POSITION - 1))

    if EXCLUDED_REGION is not None:
        for start, length in EXCLUDED_REGION:
            e_ranges.append(('EXCLUDED_REGION', start - 1, start + length - 2))
        e_ranges = remove_overlap(e_ranges)
        ranges += e_ranges

    if TARGET is not None:
        for start, length in TARGET:
            t_ranges.append(('TARGET', start - 1, start + length - 2))
        t_ranges = remove_overlap(t_ranges)
        ranges += t_ranges
        if io_range is not None:
            p_ranges += non_intersecting_parts(t_ranges, io_range)
    elif io_range is not None:
        p_ranges.append(io_range)

    if EXCLUDED_REGION is not None or TARGET is not None:
        range_gaps = remove_overlap(e_ranges + t_ranges)

    for p_range in p_ranges:
        if len(range_gaps) > 0:
            ranges += non_intersecting_parts(range_gaps, p_range)
        else:
            ranges += [p_range]

    if INCLUDED_REGION is not None:
        i_range = ''
        for start, length in INCLUDED_REGION:
            ir_start = start - 1
            ir_end = start + length - 2
            i_range = ('INCLUDED_REGION', ir_start, ir_end)
            break
        range_gaps = remove_overlap(ranges)
        ranges += non_intersecting_parts(range_gaps, i_range)
    print(ranges)
    for label, start, end in ranges:
        start_line = start // 50
        end_line = end // 50
        if start_line == end_line:
            block_coords.append((label, start, end))
        else:
            block_coords.append((label, start, (start_line + 1) * 50 - 1))
            for i in range(start_line + 1, end_line):
                block_coords.append((label, i * 50, (i + 1) * 50 - 1))
            block_coords.append((label, end_line * 50, end))
    formatted_seq = ''
    c = 0
    close = False
    block_coords = sorted(block_coords, key=lambda x: x[1])

    for i, char in enumerate(seq):
        if c < len(block_coords):
            label, block_start, block_end = block_coords[c]
            if i % 50 == 0:
                formatted_seq += str(i + 1)
                padding = 0
                while padding < 12 - len(str(i + 1)):
                    formatted_seq += "&nbsp;"
                    padding += 1
            if i == block_start:
                formatted_seq += f'<span style="color:black; background-color: {colors[label]}">'
                close = True
            if (i + 1) % 10 == 0 and (i + 1) % 50 != 0:
                formatted_seq += char
                if io_range is not None:
                    io_start = io_range[1]
                    io_end = io_range[2]
                else:
                    io_start = len(seq) + 1
                    io_end = 0
                if i == block_end and io_start < i <= io_end and label == 'TARGET':
                    formatted_seq += '</span><span style="color:black; background-color: ' + colors['INTERNAL'] + '">'
                elif i > 0 and i == ir_end and label == 'INCLUDED_REGION':
                    formatted_seq += '</span><span>'
                elif i == block_end and ir_start < i < ir_end and label != 'INCLUDED_REGION':
                    formatted_seq += '</span><span style="color:black; background-color: ' + colors['INCLUDED_REGION'] + '">'
                elif i == block_end and ir_start >= i >= ir_end:
                    formatted_seq += '</span><span>'
                formatted_seq += '&nbsp;&nbsp;'
            elif (i + 1) % 50 == 0:
                formatted_seq += f'{char}<br>\n'
            else:
                formatted_seq += char
            if i == block_end and close:
                formatted_seq += "</span>"
                c += 1
                close = False
        else:
            if i % 50 == 0:
                formatted_seq += str(i + 1)
                padding = 0
                while padding < 12 - len(str(i + 1)):
                    formatted_seq += "&nbsp;"
                    padding += 1
            if (i + 1) % 10 == 0 and (i + 1) % 50 != 0:
                formatted_seq += f'{char}&nbsp;&nbsp;'
            elif (i + 1) % 50 == 0:
                formatted_seq += f'{char}<br>\n'
            else:
                formatted_seq += char
    return f'<span style="font-family: monospace, \'Courier New\';">{formatted_seq}</span>'


def ranges_to_list(input):
    r_list = input.split()
    output = []
    for s in r_list:
        if ',' in s:
            rs = s.split(',')
        elif '-' in s:
            rs = s.split('-')
        output.append([int(x) for x in rs])
    return output


def text_monospace(text):
    html = f'<span style="font-family: monospace, \'Courier New\';">{text}</span>'
    return html


def highlight(text, color):
    html = f'<span style="color:black; background-color: {color}">{text}</span>'
    return html


def reset_values(st_values={}):
    st_values = dict((x, {}) for x in list(st_default_values.keys()))
    for key in st_values:
        st_values[key] = dict(st_default_values[key])
    return st_values


st.set_page_config(page_title="Primer3-Streamlit", page_icon=":dna:",
                   layout="centered", initial_sidebar_state="collapsed")
st.markdown("""
        <style>
               .block-container {
                    padding-top: 1rem;
                    padding-bottom: 0rem;
                }
        </style>
        """, unsafe_allow_html=True)

params = {}
st_values = reset_values()
table_th_default = 0
table_salt_default = 0
seq_args = {}
global_args = {}

for key, value in st_static_values.items():
    if key.startswith("SEQUENCE_"):
        seq_args[key] = value
    elif key.startswith("PRIMER_"):
        global_args[key] = value
p3_args = []
misprime_lib = None
misprime_lib_name = None
primer_desc = {'POSITION': 'Start',
               'LENGTH': 'Length',
               'TM': 'Tm',
               'GC_PERCENT': 'CG',
               'SELF_ANY': 'ANY',
               'SELF_END': 'SELF'
               }
primer_type = {'LEFT': 'Left Primer',
               'INTERNAL': 'Internal Oligo',
               'RIGHT': 'Right Primer'
               }
if "load_example" not in st.session_state:
    st.session_state.load_example = False
if "pick_primers" not in st.session_state:
    st.session_state.pick_primers = False
if "reset_form" not in st.session_state:
    st.session_state.reset_form = False
if "task" not in st.session_state:
    st.session_state.task = "Detection"
if "task_index" not in st.session_state:
    st.session_state.task_index = 0
if "task_help" not in st.session_state:
    st.session_state.task_help = task_help[st.session_state.task]
if "back" not in st.session_state:
    st.session_state.back = False
if "misprime_lib_index" not in st.session_state:
    st.session_state.misprime_lib_index = 0
if "table_th_index" not in st.session_state:
    st.session_state.table_th_index = table_th_default
if "table_salt_index" not in st.session_state:
    st.session_state.table_salt_index = table_th_default
primer3_main = st.empty()


def change_task():
    st.session_state.task = st.session_state["SCRIPT_TASK"]
    st.session_state.task_help = task_help[st.session_state["SCRIPT_TASK"]]
    st.session_state.task_index = st_args["SCRIPT_TASK"]["options"].index(
        st.session_state["SCRIPT_TASK"])
    return


def back_to_input():
    st.session_state.pick_primers = False
    st.session_state.back = params
    return


def reset_session_keys():
    for key in st.session_state:
        try:
            st.session_state[key] = st_default_values[key]["value"]
        except KeyError:
            continue
    st.session_state["SCRIPT_PRIMER_TM_FORMULA"] = st_args["SCRIPT_PRIMER_TM_FORMULA"]["options"][table_th_default]
    st.session_state["SCRIPT_PRIMER_SALT_CORRECTIONS"] = st_args["SCRIPT_PRIMER_SALT_CORRECTIONS"]["options"][table_salt_default]
    st.session_state["PRIMER_MISPRIMING_LIBRARY"] = st_args["PRIMER_MISPRIMING_LIBRARY"]["options"][0]
    return


with primer3_main.container():
    st.title('Primer3')
    col_1, col_2, col_3, col_4 = st.columns([3, 3, 1, 1])
    with col_1:
        st.write('Pick primers from a DNA sequence')
    with col_3:
        st.write('[About](/About)')
    with col_4:
        st.write('[Help](/Help)')
    col_1, col_2, col_3, col_4, col_5 = st.columns([4, 4, 1, 3, 3])
    with col_1:
        st.selectbox(**st_args["SCRIPT_TASK"], on_change=change_task,
                     index=st.session_state.task_index, help=st.session_state.task_help)
    with col_2:
        params["SCRIPT_SETTINGS_PRESET"] = st.selectbox(**st_args["SCRIPT_SETTINGS_PRESET"],
                                                        help="Not implemented yet", disabled=True)
    with col_4:
        st.caption('')
        st.caption('')
        if st.button('Reset form', type="secondary"):
            st.session_state.reset_form = True
    with col_5:
        st.caption('')
        st.caption('')
        if st.button('Pick Primers', key="PICK_PRIMERS", type="primary"):
            st.session_state.pick_primers = True
    tab_main, tab_g_set, tab_a_set, tab_io, tab_pw, tab_sq = st.tabs(
        ['Main', 'General Settings', 'Advanced Settings', 'Internal Oligo', 'Penalty Weights', 'Sequence Quality'])

    with tab_main:
        primer3_input = st.empty()
        with primer3_input.container():
            if st.session_state.reset_form:
                params = {}
                task = 'Detection'
                st_values = reset_values(st_values)
                reset_session_keys()
                st.session_state.table_th_index = 0
                st.session_state.table_salt_index = 0
                st.session_state.misprime_lib_index = 0
                st.session_state.load_example = False
                st.session_state.reset_form = False
            if st.session_state.back:
                for key in st_values:
                    if key in st.session_state.back and 'value' in st_values[key]:
                        st_values[key]["value"] = st.session_state.back[key]
                        params[key] = st.session_state.back[key]
            col_1, col_2 = st.columns([1, 2])
            with col_1:
                if st.session_state.task not in ["Primer_Check"]:
                    if st.button('Load Example', type="secondary"):
                        st.session_state.reset_form = False
                        st.session_state.pick_primers = False
                        st.session_state.load_example = True
                if st.session_state.load_example:
                    for key in example_values:
                        st_values[key]["value"] = example_values[key]
                params["SCRIPT_SEQUENCE_ID"] = st.text_input(
                    **st_args["SCRIPT_SEQUENCE_ID"], **st_values["SCRIPT_SEQUENCE_ID"])
                if st.session_state.task in ["Primer_Check"]:
                    params["SEQUENCE_PRIMER"] = st.text_input("Primer to test:", key="SEQUENCE_PRIMER", **st_values["SEQUENCE_PRIMER"])
            if st.session_state.task in ["Detection", "Cloning", "Sequencing", "Primer_List"]:
                with col_2:
                    st.caption('')
                    params["SCRIPT_SEQUENCE_FILE"] = st.file_uploader(
                        **st_args["SCRIPT_SEQUENCE_FILE"], **st_values["SCRIPT_SEQUENCE_FILE"], help="Not implemented yet", disabled=True)
                params["SEQUENCE_TEMPLATE"] = st.text_area(
                    **st_args["SEQUENCE_TEMPLATE"], **st_values["SEQUENCE_TEMPLATE"])
                col_1, col_2, col_3, col_4, col_5, col_6, col_7 = st.columns([6, 2, 2, 2, 3, 1, 5])
                with col_1:
                    st.markdown("Mark selected region:", help="Not implemented yet")
                with col_2:
                    if st.session_state.task in ["Detection", "Primer_List"]:
                        st.button("&lt;&nbsp;&gt;", disabled=True)
                with col_3:
                    if st.session_state.task in ["Detection", "Sequencing"]:
                        st.button("[&nbsp;]", disabled=True)
                with col_4:
                    if st.session_state.task in ["Detection", "Cloning", "Primer_List"]:
                        st.button("{&nbsp;}", disabled=True)
                with col_5:
                    st.button("Clear", disabled=True)
                with col_7:
                    st.button("Save Sequence", help="Not implemented yet", disabled=True)
                col_a, col_b, col_c = st.columns(3)
                with col_a:
                    if st.session_state.task in ["Detection", "Primer_List"]:
                        col_1, col_2, col_3 = st.columns([1, 18, 1])
                        with col_1:
                            st.caption('')
                            st.caption('')
                            st.write("&lt;&nbsp;")
                        with col_2:
                            params["SCRIPT_EXCLUDED_REGION"] = st.text_input(
                                **st_args["SCRIPT_EXCLUDED_REGION"], **st_values["SCRIPT_EXCLUDED_REGION"])
                        with col_3:
                            st.caption('')
                            st.caption('')
                            st.write("&gt;&nbsp;")
                with col_b:
                    if st.session_state.task in ["Detection", "Sequencing"]:
                        col_1, col_2, col_3 = st.columns([1, 18, 1])
                        with col_1:
                            st.caption('')
                            st.caption('')
                            st.write("[&nbsp;")
                        with col_2:
                            params["SCRIPT_TARGET"] = st.text_input(
                                **st_args["SCRIPT_TARGET"], **st_values["SCRIPT_TARGET"])
                        with col_3:
                            st.caption('')
                            st.caption('')
                            st.write("]&nbsp;")
                with col_c:
                    if st.session_state.task in ["Detection", "Cloning", "Primer_List"]:
                        col_1, col_2, col_3 = st.columns([1, 18, 1])
                        with col_1:
                            st.caption('')
                            st.caption('')
                            st.write("{&nbsp;")
                        with col_2:
                            params["SCRIPT_INCLUDED_REGION"] = st.text_input(
                                **st_args["SCRIPT_INCLUDED_REGION"], **st_values["SCRIPT_INCLUDED_REGION"])
                        with col_3:
                            st.caption('')
                            st.caption('')
                            st.write("}&nbsp;")
                if st.session_state.task in ["Detection"]:
                    col_1, col_2, col_3 = st.columns(3)
                    with col_1:
                        params["SCRIPT_DETECTION_PICK_LEFT"] = st.checkbox(
                            **st_args["SCRIPT_DETECTION_PICK_LEFT"], **st_values["SCRIPT_DETECTION_PICK_LEFT"])
                        params["SEQUENCE_PRIMER"] = st.text_input(
                            **st_args["SEQUENCE_PRIMER"], **st_values["SEQUENCE_PRIMER"])
                    with col_2:
                        params["SCRIPT_DETECTION_PICK_HYB_PROBE"] = st.checkbox(
                            **st_args["SCRIPT_DETECTION_PICK_HYB_PROBE"], **st_values["SCRIPT_DETECTION_PICK_HYB_PROBE"])
                        params["SEQUENCE_INTERNAL_OLIGO"] = st.text_input(
                            **st_args["SEQUENCE_INTERNAL_OLIGO"], **st_values["SEQUENCE_INTERNAL_OLIGO"])
                    with col_3:
                        params["SCRIPT_DETECTION_PICK_RIGHT"] = st.checkbox(
                            **st_args["SCRIPT_DETECTION_PICK_RIGHT"], **st_values["SCRIPT_DETECTION_PICK_RIGHT"])
                        params["SEQUENCE_PRIMER_REVCOMP"] = st.text_input(
                            **st_args["SEQUENCE_PRIMER_REVCOMP"], **st_values["SEQUENCE_PRIMER_REVCOMP"])

    with tab_g_set:
        params["PRIMER_PRODUCT_SIZE_RANGE"] = st.text_input(
            **st_args["PRIMER_PRODUCT_SIZE_RANGE"], **st_values["PRIMER_PRODUCT_SIZE_RANGE"])
        col_1, col_2, col_3, col_4, col_5 = st.columns([4, 7, 7, 7, 15])
        with col_1:
            st.caption('')
            st.write('[Primer Size](/Help#PRIMER_SIZE)')
        with col_2:
            params["PRIMER_MIN_SIZE"] = st.number_input(
                "Min:", min_value=1, max_value=36, value=18, step=1, key="PRIMER_MIN_SIZE")
        with col_3:
            params["PRIMER_OPT_SIZE"] = st.number_input(
                "Opt:", min_value=1, max_value=36, value=20, step=1, key="PRIMER_OPT_SIZE")
        with col_4:
            params["PRIMER_MAX_SIZE"] = st.number_input(
                "Max:", min_value=1, max_value=36, value=27, step=1, key="PRIMER_MAX_SIZE")
        with col_5:
            st.caption('')
        col_1, col_2, col_3, col_4, col_5 = st.columns([4, 7, 7, 7, 15])
        with col_1:
            st.caption('')
            st.write('[Primer TM](/Help#PRIMER_TM)')
        with col_2:
            params["PRIMER_MIN_TM"] = st.number_input(
                "Min:", min_value=0.0, max_value=100.0, value=57.0, step=0.1, key="PRIMER_MIN_TM")
        with col_3:
            params["PRIMER_OPT_TM"] = st.number_input(
                "Opt:", min_value=0.0, max_value=100.0, value=60.0, step=0.1, key="PRIMER_OPT_TM")
        with col_4:
            params["PRIMER_MAX_TM"] = st.number_input(
                "Max:", min_value=0.0, max_value=100.0, value=63.0, step=0.1, key="PRIMER_MAX_TM")
        with col_5:
            st.caption('')
        # params["PRIMER_MAX_DIFF_TM"] = st.number_input(**st_args["PRIMER_MAX_DIFF_TM"], **st_values["PRIMER_MAX_DIFF_TM"])
        # params["PRIMER_PRODUCT_SIZE_RANGE"] = st.text_input(**st_args["PRIMER_PRODUCT_SIZE_RANGE"], **st_values["PRIMER_PRODUCT_SIZE_RANGE"])
        col_1, col_2, col_3, col_4, col_5 = st.columns([4, 7, 7, 7, 15])
        with col_1:
            st.caption('')
            st.write('[Primer %GC](/Help#PRIMER_GC_PERCENT)')
        with col_2:
            params["PRIMER_MIN_GC"] = st.number_input(
                "Min:", min_value=0.0, max_value=100.0, value=20.0, step=0.1, key="PRIMER_MIN_GC")
        with col_3:
            params["PRIMER_OPT_GC"] = st.number_input(
                "Opt:", min_value=0.0, max_value=100.0, value=50.0, step=0.1, key="PRIMER_OPT_GC")
        with col_4:
            params["PRIMER_MAX_GC"] = st.number_input(
                "Max:", min_value=0.0, max_value=100.0, value=80.0, step=0.1, key="PRIMER_MAX_GC")
        with col_5:
            # params["PRIMER_GC_PERCENT"] = st.number_input(**st_args["PRIMER_GC_PERCENT"], **st_values["PRIMER_GC_PERCENT"])
            col_1, col_2, col_3 = st.columns([5, 6, 10])
            with col_1:
                st.caption('')
                st.caption('')
                st.markdown("[Fix the](/Help#SCRIPT_FIX_PRIMER_END)")
            with col_2:
                params["SCRIPT_FIX_PRIMER_END"] = st.selectbox(
                    label="_", label_visibility="hidden", options=(5, 3), key="SCRIPT_FIX_PRIMER_END")
            with col_3:
                st.caption('')
                st.caption('')
                st.markdown("[prime end of the primer](/Help#SCRIPT_FIX_PRIMER_END)",
                            **st_args["SCRIPT_FIX_PRIMER_END"], **st_values["SCRIPT_FIX_PRIMER_END"])
        col_1, col_2 = st.columns(2)
        with col_1:
            params["PRIMER_SALT_MONOVALENT"] = st.number_input(
                **st_args["PRIMER_SALT_MONOVALENT"], **st_values["PRIMER_SALT_MONOVALENT"])
            params["PRIMER_SALT_DIVALENT"] = st.number_input(
                **st_args["PRIMER_SALT_DIVALENT"], **st_values["PRIMER_SALT_DIVALENT"])
        with col_2:
            params["PRIMER_DNA_CONC"] = st.number_input(
                **st_args["PRIMER_DNA_CONC"], **st_values["PRIMER_DNA_CONC"])
            params["PRIMER_DNTP_CONC"] = st.number_input(
                **st_args["PRIMER_DNTP_CONC"], **st_values["PRIMER_DNTP_CONC"])
        st.divider()
        col_1, col_2 = st.columns(2)
        with col_1:
            params["PRIMER_MISPRIMING_LIBRARY"] = st.selectbox(
                **st_args["PRIMER_MISPRIMING_LIBRARY"], index=st.session_state.misprime_lib_index)
        st.divider()
        st.markdown("To upload or save a settings file from your local computer, choose here:",
                    help="Not implemented yet")
        col_1, col_2, col_3 = st.columns([11, 5, 4])
        with col_1:
            params["SCRIPT_SETTINGS_FILE"] = st.file_uploader(
                **st_args["SCRIPT_SETTINGS_FILE"], **st_values["SCRIPT_SETTINGS_FILE"], disabled=True)
        with col_2:
            st.caption('')
            st.caption('')
            st.caption('')
            params["Activate_Settings"] = st.button(label="Activate Settings", disabled=True)
        with col_3:
            st.caption('')
            st.caption('')
            st.caption('')
            params["Save_Settings"] = st.download_button(
                label="Save Settings", data="", disabled=True)

    with tab_a_set:
        col_1, col_2 = st.columns(2)
        with col_1:
            params["PRIMER_MAX_POLY_X"] = st.number_input(
                **st_args["PRIMER_MAX_POLY_X"], **st_values["PRIMER_MAX_POLY_X"])
            params["PRIMER_MAX_NS_ACCEPTED"] = st.number_input(
                **st_args["PRIMER_MAX_NS_ACCEPTED"], **st_values["PRIMER_MAX_NS_ACCEPTED"])
            params["PRIMER_NUM_RETURN"] = st.number_input(
                **st_args["PRIMER_NUM_RETURN"], **st_values["PRIMER_NUM_RETURN"])
            params["PRIMER_MAX_SELF_ANY"] = st.number_input(
                **st_args["PRIMER_MAX_SELF_ANY"], **st_values["PRIMER_MAX_SELF_ANY"])
        with col_2:
            params["SCRIPT_PRIMER_TM_FORMULA"] = st.selectbox(
                **st_args["SCRIPT_PRIMER_TM_FORMULA"], index=st.session_state.table_th_index)
            params["SCRIPT_PRIMER_SALT_CORRECTIONS"] = st.selectbox(
                **st_args["SCRIPT_PRIMER_SALT_CORRECTIONS"], index=st.session_state.table_salt_index)
            params["PRIMER_GC_CLAMP"] = st.number_input(
                **st_args["PRIMER_GC_CLAMP"], **st_values["PRIMER_GC_CLAMP"])
            params["PRIMER_MAX_SELF_END"] = st.number_input(
                **st_args["PRIMER_MAX_SELF_END"], **st_values["PRIMER_MAX_SELF_END"])
        col_1, col_2 = st.columns(2)
        with col_2:
            params["PRIMER_MAX_END_STABILITY"] = st.number_input(
                **st_args["PRIMER_MAX_END_STABILITY"], **st_values["PRIMER_MAX_END_STABILITY"])
        col_1, col_2 = st.columns(2)
        with col_1:
            params["PRIMER_MAX_LIBRARY_MISPRIMING"] = st.number_input(
                **st_args["PRIMER_MAX_LIBRARY_MISPRIMING"], **st_values["PRIMER_MAX_LIBRARY_MISPRIMING"])
            params["PRIMER_MAX_TEMPLATE_MISPRIMING"] = st.number_input(
                **st_args["PRIMER_MAX_TEMPLATE_MISPRIMING"], **st_values["PRIMER_MAX_TEMPLATE_MISPRIMING"])
            params["SCRIPT_PRIMER_NAME_ACRONYM_LEFT"] = st.text_input(
                **st_args["SCRIPT_PRIMER_NAME_ACRONYM_LEFT"], **st_values["SCRIPT_PRIMER_NAME_ACRONYM_LEFT"])
            params["SCRIPT_PRIMER_NAME_ACRONYM_RIGHT"] = st.text_input(
                **st_args["SCRIPT_PRIMER_NAME_ACRONYM_RIGHT"], **st_values["SCRIPT_PRIMER_NAME_ACRONYM_RIGHT"])
        with col_2:
            params["PRIMER_PAIR_MAX_LIBRARY_MISPRIMING"] = st.number_input(
                **st_args["PRIMER_PAIR_MAX_LIBRARY_MISPRIMING"], **st_values["PRIMER_PAIR_MAX_LIBRARY_MISPRIMING"])
            params["PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING"] = st.number_input(
                **st_args["PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING"], **st_values["PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING"])
            params["SCRIPT_PRIMER_NAME_ACRONYM_INTERNAL"] = st.text_input(
                **st_args["SCRIPT_PRIMER_NAME_ACRONYM_INTERNAL"], **st_values["SCRIPT_PRIMER_NAME_ACRONYM_INTERNAL"])
            params["SCRIPT_PRIMER_NAME_ACRONYM_SPACER"] = st.text_input(
                **st_args["SCRIPT_PRIMER_NAME_ACRONYM_SPACER"], **st_values["SCRIPT_PRIMER_NAME_ACRONYM_SPACER"])
        # params["PRIMER_PRODUCT_TM"] = st.number_input(**st_args["PRIMER_PRODUCT_TM"], **st_values["PRIMER_PRODUCT_TM"])
        # params["PRIMER_PRODUCT_SIZE"] = st.number_input(**st_args["PRIMER_PRODUCT_SIZE"], **st_values["PRIMER_PRODUCT_SIZE"])
        col_1, col_2, col_3 = st.columns([3, 7, 5])
        with col_1:
            params["SCRIPT_PRIMER_LIBERAL_BASE"] = st.checkbox(
                **st_args["SCRIPT_PRIMER_LIBERAL_BASE"], **st_values["SCRIPT_PRIMER_LIBERAL_BASE"])
        with col_2:
            params["SCRIPT_PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS"] = st.checkbox(
                **st_args["SCRIPT_PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS"], **st_values["SCRIPT_PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS"])
        with col_3:
            params["SCRIPT_PRIMER_LOWERCASE_MASKING"] = st.checkbox(
                **st_args["SCRIPT_PRIMER_LOWERCASE_MASKING"], **st_values["SCRIPT_PRIMER_LOWERCASE_MASKING"])
        st.divider()
        st.write("Sequencing:")
        col_1, col_2 = st.columns(2)
        with col_1:
            params["SCRIPT_SEQUENCING_LEAD"] = st.number_input(
                **st_args["SCRIPT_SEQUENCING_LEAD"], **st_values["SCRIPT_SEQUENCING_LEAD"])
            params["SCRIPT_SEQUENCING_SPACING"] = st.number_input(
                **st_args["SCRIPT_SEQUENCING_SPACING"], **st_values["SCRIPT_SEQUENCING_SPACING"])
        with col_2:
            params["SCRIPT_SEQUENCING_ACCURACY"] = st.number_input(
                **st_args["SCRIPT_SEQUENCING_ACCURACY"], **st_values["SCRIPT_SEQUENCING_ACCURACY"])
            params["SCRIPT_SEQUENCING_INTERVAL"] = st.number_input(
                **st_args["SCRIPT_SEQUENCING_INTERVAL"], **st_values["SCRIPT_SEQUENCING_INTERVAL"])
        col_1, col_2 = st.columns(2)
        with col_1:
            params["SCRIPT_SEQUENCING_REVERSE"] = st.number_input(
                **st_args["SCRIPT_SEQUENCING_REVERSE"], **st_values["SCRIPT_SEQUENCING_REVERSE"])
        st.divider()
        col_1, col_2, col_3 = st.columns(3)
        with col_1:
            st.write('Show sequence block:')
            params["SCRIPT_SEQUENCE_BLOCK_PAIRS"] = st.checkbox("Only on pair design", value=True)
            params["SCRIPT_SEQUENCE_BLOCK_FIRST"] = st.checkbox("Only on first entry", value=True)
        with col_2:
            params["SCRIPT_SHOW_INPUT"] = st.checkbox("Show original input in JSON format", value=True)
        with col_3:
            results_format = st.radio('Show output in JSON format:', options=('No', 'Flat', 'Hierarchized'), help='"Flat" is the standard Primer3 output. "Hierarchized" means outputing Primer3 results in an hierachical dictionary, which is easier to read and use in downstream applications.')
            if results_format != 'No':
                params[f'SCRIPT_SHOW_OUTPUT_{results_format.upper()}'] = True
            else:
                params['SCRIPT_SHOW_OUTPUT_FLAT'] = False
                params['SCRIPT_SHOW_OUTPUT_HIERARCHIZED'] = False
    with tab_io:
        st.caption('Not implemented yet')
    with tab_pw:
        st.caption('Not implemented yet')
    with tab_sq:
        st.caption('Not implemented yet')

global_args["PRIMER_TASK"] = primer_task[st.session_state.task]
if st.session_state.task in ["Primer_Check"]:
    params["SEQUENCE_TEMPLATE"] = params["SEQUENCE_PRIMER"]
    global_args["PRIMER_TM_FORMULA"] = 0
    global_args["PRIMER_SALT_CORRECTIONS"] = 0
if st.session_state.pick_primers and params["SEQUENCE_TEMPLATE"] != "":
    primer3_main.empty()
    params["SEQUENCE_TEMPLATE"] = params["SEQUENCE_TEMPLATE"].strip()
    illegal_chars = []
    if any(char.isdigit()for char in params["SEQUENCE_TEMPLATE"]):
        illegal_chars.append('numbers')
    if any(char.isspace() for char in params["SEQUENCE_TEMPLATE"]):
        illegal_chars.append('spaces')
    if len(illegal_chars) > 0:
        params["SEQUENCE_TEMPLATE"] = ''.join(filter(lambda char: not char.isdigit() and not char.isspace(), params["SEQUENCE_TEMPLATE"]))
        st.warning('Deleted ' + ' and '.join(illegal_chars) + ' in input sequence', icon="⚠️")
    st_values["SEQUENCE_TEMPLATE"]["value"] = params["SEQUENCE_TEMPLATE"]
    for key, value in seq_args.items():
        p3_args.append(f'{key}={value}')
    for key in params:
        if params[key] == "":
            continue
        if key.startswith("SEQUENCE_"):
            seq_args[key] = params[key]
            p3_args.append(f'{key}={params[key]}')
        if key.startswith("PRIMER_") and "_MISPRIMING_LIBRARY" not in key:
            global_args[key] = params[key]
        if params["SCRIPT_PRIMER_LIBERAL_BASE"]:
            global_args["PRIMER_LIBERAL_BASE"] = 1
        else:
            global_args["PRIMER_LIBERAL_BASE"] = 0
        if params["SCRIPT_PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS"]:
            global_args["PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS"] = 1
        else:
            global_args["PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS"] = 0
        if params["SCRIPT_PRIMER_LOWERCASE_MASKING"]:
            global_args["PRIMER_LOWERCASE_MASKING"] = 1
        else:
            global_args["PRIMER_LOWERCASE_MASKING"] = 0
    if st.session_state.task not in ["Primer_Check"]:
        regions = []
        if st.session_state.task in ["Detection", "Primer_List"]:
            regions.append("EXCLUDED_REGION")
        if st.session_state.task in ["Detection", "Sequencing"]:
            regions.append("TARGET")
        if st.session_state.task in ["Detection", "Cloning", "Primer_List"]:
            regions.append("INCLUDED_REGION")
        for r in regions:
            if params["SCRIPT_" + r] != "":
                p3_args.append("SEQUENCE_" + r + "=" + params["SCRIPT_" + r])
                seq_args["SEQUENCE_" + r] = ranges_to_list(params["SCRIPT_" + r])
        if st.session_state.task in ["Detection"]:
            if params["SCRIPT_DETECTION_PICK_LEFT"]:
                global_args["PRIMER_PICK_LEFT_PRIMER"] = 1
            else:
                global_args["PRIMER_PICK_LEFT_PRIMER"] = 0
            if params["SCRIPT_DETECTION_PICK_HYB_PROBE"]:
                global_args["PRIMER_PICK_INTERNAL_OLIGO"] = 1
            else:
                global_args["PRIMER_PICK_INTERNAL_OLIGO"] = 0
            if params["SCRIPT_DETECTION_PICK_RIGHT"]:
                global_args["PRIMER_PICK_RIGHT_PRIMER"] = 1
            else:
                global_args["PRIMER_PICK_RIGHT_PRIMER"] = 0
        global_args["PRIMER_TM_FORMULA"] = table_th[params["SCRIPT_PRIMER_TM_FORMULA"]]
        st.session_state.table_th_index = table_th[params["SCRIPT_PRIMER_TM_FORMULA"]]
        global_args["PRIMER_SALT_CORRECTIONS"] = table_salt[params["SCRIPT_PRIMER_SALT_CORRECTIONS"]]
        st.session_state.table_salt_index = table_salt[params["SCRIPT_PRIMER_SALT_CORRECTIONS"]]
        if params["PRIMER_MISPRIMING_LIBRARY"] != 'NONE':
            p3_args.append("PRIMER_MISPRIMING_LIBRARY" + "=" + params["PRIMER_MISPRIMING_LIBRARY"])
            misprime_lib = getattr(misprime_libs, params["PRIMER_MISPRIMING_LIBRARY"])
            st.session_state.misprime_lib_index = st_args["PRIMER_MISPRIMING_LIBRARY"]["options"].index(
                params["PRIMER_MISPRIMING_LIBRARY"])
    for key, value in global_args.items():
        p3_args.append(f'{key}={value}')
    p3_args = sorted(p3_args)
    p3_args.append('=')
    if 'PRIMER_PRODUCT_SIZE_RANGE' in global_args:
        global_args['PRIMER_PRODUCT_SIZE_RANGE'] = ranges_to_list(global_args['PRIMER_PRODUCT_SIZE_RANGE'])
    output_global_args = dict(global_args)
    try:
        primer3_results = design_primers(
            seq_args, global_args, misprime_lib=misprime_lib)
        primers = hierarchize(primer3_results)
    except Exception as e:
        st.exception(e)
        st.subheader('Original input')
        st.json({'seq_args': seq_args, 'global_args': output_global_args, 'misprime_lib': misprime_lib})
        try:
            st.subheader('Original output')
            st.write(primer3_results)
        except Exception:
            pass

    sequence_block_params = {'seq': params["SEQUENCE_TEMPLATE"]}
    for region in ['SEQUENCE_EXCLUDED_REGION', 'SEQUENCE_TARGET', 'SEQUENCE_INCLUDED_REGION']:
        bp = region.replace('SEQUENCE_', '')
        if region in seq_args:
            sequence_block_params[bp] = seq_args[region]
    if 'WARNING' in primers:
        st.warning(primers['WARNING'], icon="⚠️")
    st.title('Primer3 Results')
    st.button('< Back', on_click=back_to_input, type="secondary")
    if len(primers['PRIMERS']) == 0:
        st.warning("No Primers found", icon="⚠️")
    for p, primer in enumerate(primers['PRIMERS']):
        primer_number = params['SCRIPT_PRIMER_NAME_ACRONYM_SPACER']
        if p > 0:
            primer_number = f'{primer_number}{p+1}{primer_number}'
        if 'PAIR' in primer:
            st.subheader(f'Pair: {p+1}')
        for pt in primer_type:
            if pt in primer:
                sequence_block_params[f'{pt}_POSITION'] = primer[pt]['POSITION']
                sequence_block_params[f'{pt}_LENGTH'] = primer[pt]['LENGTH']
                col_1, col_2 = st.columns(2)
                with col_1:
                    if params["SCRIPT_SEQUENCE_ID"] != '':
                        st.write(highlight(
                            f'{primer_type[pt]} {p+1}: {params["SCRIPT_SEQUENCE_ID"]}{primer_number}' + params[f'SCRIPT_PRIMER_NAME_ACRONYM_{pt}'], colors[pt]), unsafe_allow_html=True)
                    else:
                        st.write(highlight(
                            f'{primer_type[pt]} {p+1}', colors[pt]), unsafe_allow_html=True)
                with col_2:
                    st.write('Sequence: ' + text_monospace(primer[pt]['SEQUENCE']), unsafe_allow_html=True)
                cols = st.columns(len(primer_desc.keys()))
                for c, key in enumerate(primer_desc):
                    with cols[c]:
                        if key in ['TM', 'GC_PERCENT']:
                            st.write(f"{primer_desc[key]}: {primer[pt][key]:.1f}")
                        else:
                            st.write(f"{primer_desc[key]}: {primer[pt][key]}")
        if 'PAIR' in primer:
            col_1, col_2, col_3 = st.columns(3)
            with col_1:
                st.write(f"Product Size: {primer['PAIR']['PRODUCT_SIZE']}")
            with col_2:
                st.write(f"Pair Any: {primer['PAIR']['COMPL_ANY']}")
            with col_3:
                st.write(f"Pair End: {primer['PAIR']['COMPL_END']}")
        if 'SEQUENCE_EXCLUDED_REGION' in seq_args:
            reg_num = len(seq_args['SEQUENCE_EXCLUDED_REGION'])
            for i, col in enumerate(st.columns(reg_num)):
                if reg_num == 1:
                    n = ''
                else:
                    n = f' {i + 1}'
                rg, le = seq_args['SEQUENCE_EXCLUDED_REGION'][i]
                with col:
                    st.write(highlight(f'Excluded region{n}:', '#FB6A6A') + f'&nbsp;&nbsp;Start: {rg}; Length: {le}', unsafe_allow_html=True)
        if 'SEQUENCE_TARGET' in seq_args:
            reg_num = len(seq_args['SEQUENCE_TARGET'])
            for i, col in enumerate(st.columns(reg_num)):
                if reg_num == 1:
                    n = ''
                else:
                    n = f' {i + 1}'
                rg, le = seq_args['SEQUENCE_TARGET'][i]
                with col:
                    st.write(highlight(f'Target{n}:', '#0DF20D') + f'&nbsp;&nbsp;Start: {rg}; Length: {le}', unsafe_allow_html=True)
        if 'SEQUENCE_INCLUDED_REGION' in seq_args:
            reg_num = len(seq_args['SEQUENCE_INCLUDED_REGION'])
            for i, col in enumerate(st.columns(reg_num)):
                if reg_num == 1:
                    n = ''
                else:
                    n = f' {i + 1}'
                rg, le = seq_args['SEQUENCE_INCLUDED_REGION'][i]
                with col:
                    st.write(highlight(f'Included region{n}:', '#83FCFC') + f'&nbsp;&nbsp;Start: {rg}; Length: {le}', unsafe_allow_html=True)
        if p == 0 and 'PAIR' in primer:
            print(sequence_block_params)
            st.write(sequence_block(**sequence_block_params), unsafe_allow_html=True)
        elif p > 0 and not params["SCRIPT_SEQUENCE_BLOCK_FIRST"]:
            st.write(sequence_block(**sequence_block_params), unsafe_allow_html=True)
        if 'PAIR' not in primer and not params["SCRIPT_SEQUENCE_BLOCK_PAIRS"]:
            if p == 0:
                st.write(sequence_block(**sequence_block_params), unsafe_allow_html=True)
            elif p > 0 and not params["SCRIPT_SEQUENCE_BLOCK_FIRST"]:
                st.write(sequence_block(**sequence_block_params), unsafe_allow_html=True)
        st.divider()
    st.subheader('Statistics')
    for primer_type, explain in primers['EXPLAIN'].items():
        if primer_type == 'PRIMER_LEFT':
            st.write(f"Left Primer: {explain}")
        if primer_type == 'PRIMER_INTERNAL':
            st.write(f"Internal Oligo: {explain}")
        if primer_type == 'PRIMER_RIGHT':
            st.write(f"Right Primer: {explain}")
        if primer_type == 'PRIMER_PAIR':
            st.write(f"Primer Pairs: {explain}")
    st.divider()
    if params["SCRIPT_SHOW_INPUT"]:
        st.subheader('Original input')
        ml = None
        if params["PRIMER_MISPRIMING_LIBRARY"] != 'NONE':
            ml = params["PRIMER_MISPRIMING_LIBRARY"]
        with st.expander('For `primer3-py`', expanded=False):
            st.json({'seq_args': seq_args, 'global_args': output_global_args, 'misprime_lib': ml})
        with st.expander('For `primer3`', expanded=False):
            st.code('\n'.join(p3_args), language=None)
    if params["SCRIPT_SHOW_OUTPUT_FLAT"]:
        st.subheader('Original output')
        with st.expander('JSON Flat', expanded=False):
            st.json(primer3_results)
    if params["SCRIPT_SHOW_OUTPUT_HIERARCHIZED"]:
        st.subheader('Original output')
        with st.expander('JSON Hierarchized', expanded=False):
            st.json(primers)
