package SimulationEnv;

import DataProvider.DataReader;
import SensorNetwork.*;
import PlacementMethods.Optimization;
import PlacementMethods.SimulationMode;

import java.io.IOException;

import SimulationLevel.*;


public class PackageSimRun {

	public static class Parameters {
		final static String network = "C:\\ARRS\\input\\250.csv";
		final static String nodeattributes =  "C:\\ARRS\\input\\node_attributes250.csv";
	}


	public static void main(String[] args) throws IOException {

		SensorNetwork sensor_Sensor_network = DataReader.readEdges(Parameters.network);
		sensor_Sensor_network = DataReader.readNodes(Parameters.nodeattributes, sensor_Sensor_network);
		Optimization.Optimization(sensor_Sensor_network, 50, 1000, PackageSimLevel.LWRL, SimulationMode.EQUALITY);
	}
}