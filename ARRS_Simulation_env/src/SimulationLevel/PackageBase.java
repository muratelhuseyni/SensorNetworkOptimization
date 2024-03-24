package SimulationLevel;

import SensorNetwork.*;

public abstract class PackageBase {
	protected SensorNode origin;
	protected SensorNode actnode;
	
	public abstract void PackageSim();
	
	public PackageBase(SensorNode origin) {
		this.origin=origin;
		this.actnode=origin;
	}

	public SensorNode getSensorName() {
		return actnode;
	}
	
}
