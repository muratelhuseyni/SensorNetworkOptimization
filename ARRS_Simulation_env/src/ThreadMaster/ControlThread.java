package ThreadMaster;

import SensorNetwork.SensorType;
import SimulationLevel.PackageBase;

public abstract class ControlThread extends Thread {
	private Thread workerThread;
	private String threadId;
	
	public ControlThread(String id) {
		this.threadId = id;
	}

	public static class ParallelPackages extends ControlThread {

		private PackageBase message;

		public ParallelPackages(String id, PackageBase message) {
			super(id);
			this.message = message;
		}

		@Override
		public void run() {
			while (!message.getSensorName().getNodeFormat().equals(SensorType.GATEWAY)) {
				message.PackageSim();
			}
		}
	}
}
