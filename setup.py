import setuptools

setuptools.setup(
    name="SpaGFT",
    version="0.0.1",
    author="Jixin Liu",
    packages=['SpaGFT'],
    author_email="frankliu210@163.com",
    description='''SpaGFT is a python package to analyze spatial transcriptomcs
                 data''',
    url="https://github.com/jxLiu-bio/SpaGFT",
    project_urls={
        "Bug Tracker": "https://github.com/jxLiu-bio/SpaGFT/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
