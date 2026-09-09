# nxsa_AI
Clustering try with XMM Newton event files.

## Preparing the python environment

* Install or compile the last Python 3.11:\
Official page: https://www.python.org

* Download the project:\
Your current folder will be called \<ROOT\>\
git clone https://github.com/davidglt/nxsa_AI.git

* Download test control dataset files: (for Linux or git bash, you can use this script)\
cd \<ROOT\>/nxsa_AI/utils\
./download_dataset_files.sh

* Create a new Python virtual environment for this project:\
cd \<ROOT\>/nxsa_AI\
<PATH_TO_YOUR_PYTHON_BIN>/python3 -m venv .venv

* Activate the environment:\
cd \<ROOT\>/nxsa_AI\
Linux: source .venv/bin/activate
Windows: .venv/scripts/activate.bat

* Import pip requirements:\
cd \<ROOT\>/nxsa_AI\
pip install -r requirements/requirements.txt

* Run jupyter-notebook and enjoy it!\
cd \<ROOT\>/nxsa_AI/src\
jupyter-notebook

## Authors

* **David Gonzalez** - [nxsa_AI](https://github.com/nxsa_AI)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the GPLv3 License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thanks to all XMM people in special to:\
**Jose Vicente Perea**\
**Pedro Rodriguez**\
**Maria Santos**\
**Norbert Schartel**
