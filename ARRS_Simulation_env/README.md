# Sensor Network Simulation

This project encompasses a comprehensive simulation environment for sensor networks, focusing on package transmission strategies, network optimization, and sensor node management. It leverages Java to simulate various package transmission protocols within a sensor network, evaluate optimization strategies for sensor placement, and provide insights into network behavior under different simulation modes.

## Features

- **Package Simulation**: Simulates different package transmission protocols across the sensor network.
- **Network Optimization**: Implements optimization strategies to determine optimal placement of gateways.
- **Data Reading**: Facilitates network configuration through CSV files, enabling easy setup and customization of the sensor network.
- **Concurrency Management**: Utilizes multithreading to simulate concurrent package transmission across the network.

## Getting Started

### Prerequisites

- Java 11 or newer

### Running the Simulation

Execute the main simulation runner class:

```sh
java -cp target/classes SimulationEnv.PackageSimRun
```

## Usage

To configure your simulation, edit the CSV files in the `input` directory according to your network setup:

- `250.csv`: Defines the sensor nodes and their connections.
- `node_attributes250.csv`: Specifies attributes for each sensor node, such as type and location.

After configuration, run the simulation as described above to evaluate the sensor network's performance under various conditions and optimization strategies.
