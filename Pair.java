
public class Pair {
	public int x;
	public int y;
	
	public Pair(int x, int y) {
		this.x = x;
		this.y = y;
	}
	
	
	/*
	 * Used this source for equals and hashcode override https://stackoverflow.com/questions/28953856/accessing-a-hashmap-value-with-a-custom-pair-object-as-a-key
	 */
	@Override
	public boolean equals(final Object o) {
		if (this == o) {
			return true;
		}
		if (!(o instanceof Pair)) {
			return false;
		}
		
		final Pair pair = (Pair) o;
		if (x != pair.x) {
			return false;
		}
		if (y != pair.y) {
			return false;
		}

		return true;
	}
	
	@Override
	public int hashCode() {
		int result = x;
		result = 31 * result + y;
		return result;
	}
}
