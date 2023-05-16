[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://primer3.streamlit.app/)

# Primer3-Streamlit

**Primer3-Streamlit** is a web app interface for Primer3 that relies on [primer3-py](https://github.com/libnano/primer3-py) and [Streamlit](https://github.com/streamlit/streamlit), and is modeled after [Primer3Plus](https://github.com/primer3-org/primer3plus).

The immediate goal is to achieve feature parity with Primer3Plus 1.1.0. Features introduced on later versions will be implemented in the near future.

Preset loading, region symbols in sequence, FASTA support and upload/download features are not implemented yet. Task loading is partially implemented. Otherwise, this web app is fully functional for primer design.


## Installation

```bash
# Download the repository from github
git clone https://github.com/pessoa-am/Primer3-StreamLit.git
cd Primer3-Streamlit
# Install dependencies
pip install -r requirements.txt
# Run with StreamLit
streamlit run Primer3-Streamlit.py
```

## Copying or Reusing

This project is free software licensed under the terms of the [GNU Affero General Public License, version 3](https://www.gnu.org/licenses/agpl-3.0.txt).

 `primer3_st_args.py`, `primer3_st_help.py` and `misprime_libs.py` are based on Primer3Plus and thus licensed under the terms of under the terms of the [GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.txt).
