import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import com.jme3.math.Vector3f;

public class Vertex {
	Vector3f pos;
	List<Vertex> neighbors = new ArrayList<Vertex>();
	double f;
	double g = 0d;
	Vector3f gradF = Vector3f.ZERO;
	Vector3f gradG = Vector3f.ZERO;
	boolean finalized = false;
	State state;
	
	@Override
	public boolean equals(Object o) {
		if(this == o) return true;
		if(o == null || getClass() != o.getClass()) return false;
		Vertex vertex = (Vertex) o;
		return pos.equals(vertex.pos);
	}
	
	@Override
	public int hashCode() {
		return Objects.hashCode(pos);
	}
}
