
## Getting Started with QuantymDotDesigner

### Installing with pip (recommended)

To install the package with pip, simply run the following command:

```bash
pip install git+https://github.com/valentinjohn/quantum-dot-designer@main
```

### Installing from a local source

Alternatively, you can clone the repository and install the package from the local source:


###### 1. Clone the Repository

To clone the repository to your local machine, open a terminal and run the following command:

```bash
git clone https://github.com/valentinjohn/quantum-dot-designer.git
```

###### 2. Set Up a New Conda Environment

To ensure that you have all the required dependencies and to avoid any conflicts, it's a good idea to set up a new Conda environment specifically for this project.

First, make sure you have [Anaconda](https://www.anaconda.com/products/distribution) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed.

Create a new Conda environment:

```bash
conda create --name dot_designer python=3.11
```

Activate the newly created environment:

```bash
conda activate dot_designer
```

###### 3. Install Dependencies

First navigate to the cloned repository
```bash
cd my_local_path
```
where you have to replace my_local_path with the path of the cloned repository.
 
Install the required dependencies from the `requirements.txt` file:

```bash
pip install -r requirements.txt
```

You're all set! Now you can start working on the project within your newly created Conda environment.
