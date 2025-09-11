# neuron_morphology_tools

Algorithms and analysis tools for networkx graph objects representing individual neurons output through the mesh processing pipeline outline in the following study: https://www.biorxiv.org/content/10.1101/2023.03.14.532674v1

## 🚀 Getting Started with Docker

You have two options for setting up the Docker environment:

### 1. Pull the Prebuilt Image

```bash
docker pull celiib/neuron_morphology:latest
```

### 2. Build the Image Locally

Navigate to the `/docker` folder and build the image using `docker-compose`:

```bash
cd docker
docker-compose build
```

---

## ▶️ Running the Environment

1. **Create a `./notebooks` directory** in the repository root (this will be mounted into the container for your work):

   ```bash
   cd docker
   mkdir notebooks
   ```

2. **Start the container** from the `/docker` folder:

   ```bash
   docker-compose up
   ```

   This will launch **JupyterLab** at [http://localhost:8890](http://localhost:8890).

---

## ⚙️ Installing the Repository Inside the Container

Once inside JupyterLab:

1. Open a **terminal** in the JupyterLab interface.
2. Run the following command to install the repository in editable mode:
   ```bash
   pip3 install -e ./neuron_morphology_tools
   ```

---

## 📓 Running Example Notebooks

After installation, you can open and run any of the example notebooks in the `Applications` folder to explore the functionality of the tools.

---

## 🧰 Requirements

- [Docker](https://docs.docker.com/get-docker/)
- [Docker Compose](https://docs.docker.com/compose/)

---

## 📝 Notes

- The `./notebooks` directory is bind-mounted so your work will persist outside the container.
- The container exposes JupyterLab on port **8890**.

---
