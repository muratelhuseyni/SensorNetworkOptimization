package SimulationLevel;

import java.util.HashSet;
import java.util.Set;

public enum PackageSimLevel {
	LWRL, MLRW, RW, IWRW;

	public static class PackageCounter extends PackageBase {

		private Set<String> messageOrigins = new HashSet<>();
		private PackageBase originalMessage;

		public PackageCounter(PackageBase message) {
			super(message.getSensorName());
			this.originalMessage = message;
			if (message.getSensorName() != null && message.getSensorName().getId() != null) {
				messageOrigins.add(message.getSensorName().getId());
			}
		}

		@Override
		public void PackageSim() {
			originalMessage.PackageSim();
			setActNode();
		}

		private void setActNode() {
			this.actnode = originalMessage.getSensorName();
		}

	}
}
