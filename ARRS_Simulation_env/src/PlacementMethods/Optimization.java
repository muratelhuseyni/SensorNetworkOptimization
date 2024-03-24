package PlacementMethods;

import PackageSimulator.PackageSimulator;
import SimulationLevel.PackageSimLevel;
import SensorNetwork.*;

public class Optimization {

	public static SensorNetwork Optimization(SensorNetwork net, int target, int accuracy, PackageSimLevel protocol, SimulationMode obj) {
		net.setAllNodeToSensor();
		while (net.getGateways().size() < target) {
			SensorNode actual_best_gateway = null;
			long all_time_best_solution = Long.MAX_VALUE;
			for (SensorNode n : net.getSensors()) {
				resetLoadForAllNodes(net);
				n.setType(SensorType.GATEWAY);
				PackageSimulator.RunMessages.RWIteration(net, accuracy, protocol);

				long actual = calculateObjective(net, obj);
				if (all_time_best_solution > actual) {
					all_time_best_solution = actual;
					actual_best_gateway = n;
				}
				n.setType(SensorType.SENSOR);
			}
			if (actual_best_gateway != null) {
				actual_best_gateway.setType(SensorType.GATEWAY);
				net.removeGatewayFromSensorList(actual_best_gateway);
				net.addGateway(actual_best_gateway);
				System.out.println("Actual best solution size: " + net.getGateways().size());
				System.out.println("Next Gateway: " + actual_best_gateway.getId() + " actual objective value: " + all_time_best_solution);
			}
		}
		return net;
	}

	private static void resetLoadForAllNodes(SensorNetwork net) {
		for (SensorNode i : net.getNodes()) {
			i.setLoad(0);
		}
	}

	private static long calculateObjective(SensorNetwork net, SimulationMode obj) {
		switch (obj) {
			case BATTERYSAVING:
				return net.getMaxLoad();
			case EQUALITY:
				return net.getBalance();
			default:
				return Long.MAX_VALUE;
		}
	}
}
