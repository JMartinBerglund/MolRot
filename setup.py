import setuptools

# Load the long_description from README.md
with open("README.md", "r") as fh:
    long_description = fh.read()

    setuptools.setup(
        name="MolRot",
        version="0.0.1",
        author="Jan Martin Berglund",
        author_email="j.martin.berglund@gmail.com",
        description="A package that performs tomography and interferometry on the rotational states of diatomic molecular ions, using off-resonance fs laser pulses.",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="",
        packages=setuptools.find_packages(include=['MolRot', 'Tutorial']),
        install_requires=[],
        setup_requires=['pytest-runner'],
        tests_require=['pytest==4.4.1'],
        test_suite='tests',
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ],
        python_requires='>=3.6',
    )
