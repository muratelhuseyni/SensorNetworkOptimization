package PackageSimulator;

import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import SensorNetwork.*;
import SimulationLevel.PackageSimLevel;
import SimulationLevel.Packages;
import ThreadMaster.ControlThread;
import ThreadMaster.IterationBase;

public class PackageSimulator {
	public static double Simulation(SensorNetwork net, LinkedList<SensorNode> nodeset, int samplesize) {
		for(int j=0;j<samplesize;j++) {

			net = GenerateG(net); //Generate G'

			for(SensorNode n:nodeset) {
				DFS(n,net);
			}
			net.empty();
		}
		net.getNodes().forEach(n -> n.finalizefv(samplesize));

		return net.getFvSum();

	}

	public static SensorNetwork GenerateG(SensorNetwork net) {
		Iterator<SensorConnection> it=net.getEdges().iterator();
		while(it.hasNext()) {
			SensorConnection e=(SensorConnection)it.next();
			if(Math.random()<e.getWeight()) {
				e.setAlive(true);
			}else {
				e.setAlive(false);
			}
		}
		return net;
	}


	public static void DFS(SensorNode v, SensorNetwork net) {
		Stack<SensorNode> s=new Stack<SensorNode>();
		s.push(v);
		while(!s.isEmpty()){
			SensorNode act=s.pop();
			if(act.getState()!= SensorInfo.VISITED) {
				act.setState(SensorInfo.VISITED);
				act.addFv(1);

				for(SensorConnection e: act.getOutlist()) {
					if(e.isAlive()) {
						s.push(e.getIn());
					}
				}

				for(SensorConnection e: act.getInlist()) {
					if(e.isAlive()) {
						s.push(e.getOut());
					}
				}
			}

		}
	}

	public static class RunMessages {
		public static void RunMessageRandomWalk(SensorNetwork net, int iterationnumber, PackageSimLevel mtype) {
			ExecutorService executor = Executors.newFixedThreadPool(10);
			for(int i=0;i<iterationnumber;i++) {
				for(SensorNode n: net.getSensors()) {
					String id = n.getId()+":"+Integer.toString(i);
					switch(mtype) {
						case RW: executor.execute( new ControlThread.ParallelPackages(id, new Packages.RWPackage(n))); break;
						case MLRW:executor.execute( new ControlThread.ParallelPackages(id, new Packages.MLRWPackage(n))); break;
						case LWRL:executor.execute( new ControlThread.ParallelPackages(id, new Packages.LWRWPackage(n))); break;
						case IWRW: executor.execute(new ControlThread.ParallelPackages(id, new Packages.WRWPackage(n)));
						default: break;
					}

				}
			}
			executor.shutdown();
		}

		public static void RWIteration(SensorNetwork net, int measurements, PackageSimLevel mtype) {
				new IterationBase("0",net, measurements).start(mtype);
		}
	}
}
