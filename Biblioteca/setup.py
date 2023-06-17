from setuptools import setup, find_packages

with open("README.md", "r") as rm:
    long_description = rm.read()

setup(
    name = "VLCDA",
    version = "1.0.4",
    license='MIT',
    description = "Modulo Python dedicado a evolução de circuitos lógicos através de algoritmos genéticos",
    author = "Humberto Távora, Stefan Blawid, Yasmin Maria",
    author_email = 'humbertocctg@gmail.com',
    long_description=long_description,
    long_description_content_type = "text/markdown",
    url = 'https://github.com/HumbertoTavora/Logic-Circuits-Evolution',
    download_url = 'https://github.com/HumbertoTavora/Logic-Circuits-Evolution',
    keywords = ['Genetic', 'Algorithm', 'Logic Circuits', 'Evolution', 'Optimization'],
    classifiers=['Programming Language :: Python :: 3','Intended Audience :: Developers'],
    packages=find_packages()
    )

