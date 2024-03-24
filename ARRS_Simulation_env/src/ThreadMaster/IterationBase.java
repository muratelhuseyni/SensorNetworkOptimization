package ThreadMaster;

import SensorNetwork.SensorNode;
import SensorNetwork.SensorType;

import java.util.ArrayList;
import java.util.List;

import SimulationLevel.*;
import SensorNetwork.SensorNetwork;

public class IterationBase extends ControlThread {

	private SensorNetwork sensorNetwork;
	private int measurements;

	public IterationBase(String id, SensorNetwork sensorNetwork, int measurements) {
		super(id);
		this.sensorNetwork = sensorNetwork;
		this.measurements = measurements;
	}

	public void start(PackageSimLevel mtype) {
		try {
			for (int i = 0; i < measurements; i++) {
				List<PackageBase> messages = new ArrayList<>();
				sensorNetwork.getSensors().forEach(sensor -> messages.add(createMessage(mtype, sensor)));
				Thread t = new Thread(new IterationWorker(messages), "iter_" + i);
				t.start();
				t.join();
			}
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			System.out.println("Interrupted during execution");
			e.printStackTrace();
		} catch (Exception e) {
			System.out.println("An error occurred");
			e.printStackTrace();
		}
	}

	private PackageBase createMessage(PackageSimLevel mtype, SensorNode sensor) {
		return switch (mtype) {
			case RW -> new Packages.RWPackage(sensor);
			case MLRW -> new Packages.MLRWPackage(sensor);
			case LWRL -> new Packages.LWRWPackage(sensor);
			case IWRW -> new Packages.WRWPackage(sensor);
			default -> throw new IllegalArgumentException("Unsupported MessageType: " + mtype);
		};
	}

	class IterationWorker extends Thread {

		private List<PackageBase> messages;

		public IterationWorker(List<PackageBase> messages) {
			this.messages = messages;
		}

		@Override
		public void run() {
			while (!messages.isEmpty()) {
				messages.removeIf(message -> {
					message.PackageSim();
					return message.getSensorName().getNodeFormat().equals(SensorType.GATEWAY);
				});
			}
		}
	}
}
