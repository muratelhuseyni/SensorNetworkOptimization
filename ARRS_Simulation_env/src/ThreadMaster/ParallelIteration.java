package ThreadMaster;

import SensorNetwork.SensorType;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import SimulationLevel.PackageBase;
import SimulationLevel.Packages;
import SensorNetwork.SensorNetwork;


import java.util.Iterator;

public class ParallelIteration extends ControlThread {

	private SensorNetwork net_data;
	private int result_data;

	public ParallelIteration(String id, SensorNetwork net_data, int result_data) {
		super(id);
		this.net_data = net_data;
		this.result_data = result_data;
	}

	@Override
	public void run() {
		System.out.println("Starting main thread");
		for (int i = 0; i < result_data; i++) {
			try {
				processIteration(i);
			} catch (InterruptedException e) {
				Thread.currentThread().interrupt();
				System.out.println("Thread was interrupted, failing to complete operation");
				return;
			}
		}
	}

	private void processIteration(int iteration) throws InterruptedException {
		List<PackageBase> collected_data = net_data.getSensors().stream()
				.map(Packages.RWPackage::new)
				.collect(Collectors.toCollection(ArrayList::new));

		Thread workerThread = new Thread(new MessageWorker(collected_data), "iter_" + iteration);
		workerThread.start();
		workerThread.join();
	}
}


class MessageWorker extends Thread {

	private List<PackageBase> messages;

	public MessageWorker(List<PackageBase> messages) {
		this.messages = messages;
	}

	@Override
	public void run() {
		System.out.println("Starting worker " + this.getName());
		// Use an Iterator to safely remove items from the list while iterating
		Iterator<PackageBase> iterator = messages.iterator();
		while (iterator.hasNext()) {
			PackageBase message = iterator.next();
			message.PackageSim();
			if (message.getSensorName().getNodeFormat().equals(SensorType.GATEWAY)) {
				iterator.remove();
			}
		}
		System.out.println("Ending worker " + this.getName());
	}
}

	
