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
from primer3_st_help import st_help


st.set_page_config(page_title="Primer3 - Help", page_icon=":dna:",
                   layout="wide", initial_sidebar_state="collapsed")

st.markdown("""
        <style>
               .block-container {
                    padding-top: 1rem;
                    padding-bottom: 0rem;
                }
        </style>
        """, unsafe_allow_html=True)

st.title("Primer3 - Help")

for entry in st_help:
    for widget, args in entry.items():
        getattr(st, widget)(**args)
