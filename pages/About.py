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


import sys
sys.path.append('..')

import streamlit as st
from primer3_st_help import st_about


st.set_page_config(page_title="Primer3 - About", page_icon=":dna:",
                   layout="wide", initial_sidebar_state="collapsed")

st.markdown("""
        <style>
               .block-container {
                    padding-top: 1rem;
                    padding-bottom: 0rem;
                }
        </style>
        """, unsafe_allow_html=True)

st.title("Primer3 - About")

st.header("Primer3-Streamlit")

st.markdown("**Primer3-Streamlit** is a web app interface for Primer3 that relies on [primer3-py](https://github.com/libnano/primer3-py) and [Streamlit](https://github.com/streamlit/streamlit), and is modeled after [Primer3Plus](https://github.com/primer3-org/primer3plus).\n\nThe immediate goal is to achieve feature parity with Primer3Plus 1.1.0. Features introduced on later versions will be implemented in the near future.\n\nPreset loading, region symbols in sequence, FASTA support and upload/download features are not implemented yet. Task loading is partially implemented. Otherwise, this web app is fully functional for primer design.")

st.subheader("Copying or Reusing")

st.markdown("Primer3-Streamlit is free software licensed under the terms of the [GNU Affero General Public License, version 3](https://www.gnu.org/licenses/agpl-3.0.txt).\n\n `primer3_st_args.py`, `primer3_st_help.py` and `misprime_libs.py` are based on Primer3Plus and thus licensed under the terms of under the terms of the [GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.txt). Original about page from Primer3Plus 1.1.0 bellow.")

for entry in st_about:
    for widget, args in entry.items():
        getattr(st, widget)(**args)
