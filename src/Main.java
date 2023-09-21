import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.TreeMap;
import java.util.stream.Collectors;

import com.jme3.app.SimpleApplication;
import com.jme3.light.AmbientLight;
import com.jme3.light.DirectionalLight;
import com.jme3.material.Material;
import com.jme3.math.ColorRGBA;
import com.jme3.math.Triangle;
import com.jme3.math.Vector3f;
import com.jme3.scene.Geometry;
import com.jme3.scene.Mesh;
import com.jme3.scene.Spatial;
import com.jme3.scene.debug.Arrow;
import com.jme3.scene.shape.Sphere;
import com.jme3.system.AppSettings;

public class Main extends SimpleApplication{
	
	public static final double THRESHOLD = 1e-4;
	public static final double GRID_WIDTH = 0.15d;

	public static void main(String[] args) {
		Main main = new Main();
		main.start();
	}
	
	@Override
	public void simpleInitApp() {
		Spatial mesh = assetManager.loadModel("sphere.obj");
		Material mat = new Material(assetManager, "Common/MatDefs/Light/Lighting.j3md");
        mat.setColor("Diffuse", ColorRGBA.Gray);
        mat.setColor("Ambient", ColorRGBA.White.mult(0.1f));
        mat.setBoolean("UseMaterialColors", true);
        mesh.setMaterial(mat);
        rootNode.attachChild(mesh);
        
        List<Vertex> verticesList = new ArrayList<Vertex>();
        
        if(mesh instanceof Geometry) {
        	Geometry geo = (Geometry) mesh;
        	Mesh m = geo.getMesh();
        	
        	for(int i = 0; i < m.getTriangleCount(); i++) {
        		Triangle t = new Triangle();
        		m.getTriangle(i, t);
        		
        		Vertex v1 = findOrCreateVertex(t.get1(), verticesList);
        		Vertex v2 = findOrCreateVertex(t.get2(), verticesList);
        		Vertex v3 = findOrCreateVertex(t.get3(), verticesList);
        		
        		if(!v1.neighbors.contains(v2)) v1.neighbors.add(v2);
        		if(!v1.neighbors.contains(v3)) v1.neighbors.add(v3);
        		if(!v2.neighbors.contains(v1)) v2.neighbors.add(v1);
        		if(!v2.neighbors.contains(v3)) v2.neighbors.add(v3);
        		if(!v3.neighbors.contains(v1)) v3.neighbors.add(v1);
        		if(!v3.neighbors.contains(v2)) v3.neighbors.add(v2);
        	}
        	
        	computeGeodesicDistance(verticesList, verticesList.get(0));
        	
        	TreeMap<Double, List<Vertex>> rows = groupByRows(verticesList);
        	
        	computeG(verticesList, 1000);
        	
        	for(Map.Entry<Double, List<Vertex>> mapEntry : rows.entrySet()) {
        		List<Vertex> list = mapEntry.getValue();
        		list.sort((v1, v2) -> Double.compare(v1.g, v2.g));
        	}
        	
//        	resampleVertices(rows);
        	
        	int index = 0;
        	for(Map.Entry<Double, List<Vertex>> entry : rows.entrySet()) {
        		List<Vertex> r = entry.getValue();
        		
        		System.out.println("row size: " + r.size());
        		for(int i = 0; i < r.size(); i++) {
        			System.out.println("vertex: " + r.get(i).pos.toString() + " f: " + r.get(i).f + " gradient of g: " + r.get(i).gradG + " gradient of f: " + r.get(i).gradF);
        			Vector3f arrowVec = new Vector3f(r.get(i).gradF).scaleAdd(0.25f, Vector3f.ZERO);
        			Arrow arrow = new Arrow(arrowVec);
        			Geometry gArrow = new Geometry("f vector", arrow);
        			Material gmat = new Material(assetManager, "Common/MatDefs/Misc/Unshaded.j3md");
        		    gmat.getAdditionalRenderState().setWireframe(true);
        		    gmat.getAdditionalRenderState().setLineWidth(4);
        		    gmat.setColor("Color", ColorRGBA.Green);
        		    gArrow.setMaterial(gmat);
        		    rootNode.attachChild(gArrow);
        		    gArrow.setLocalTranslation(r.get(i).pos);
        		    
        		    arrowVec = new Vector3f(r.get(i).rotatedGradF).scaleAdd(0.25f, Vector3f.ZERO);
        			arrow = new Arrow(arrowVec);
        			gArrow = new Geometry("rotated f vector", arrow);
        			gmat = new Material(assetManager, "Common/MatDefs/Misc/Unshaded.j3md");
        		    gmat.getAdditionalRenderState().setWireframe(true);
        		    gmat.getAdditionalRenderState().setLineWidth(4);
        		    gmat.setColor("Color", ColorRGBA.Red);
        		    gArrow.setMaterial(gmat);
        		    rootNode.attachChild(gArrow);
        		    gArrow.setLocalTranslation(r.get(i).pos);
        		    
//        		    arrowVec = new Vector3f(r.get(i).rotatedGradF2).scaleAdd(0.25f, Vector3f.ZERO);
//        			arrow = new Arrow(arrowVec);
//        			gArrow = new Geometry("g vector", arrow);
//        			gmat = new Material(assetManager, "Common/MatDefs/Misc/Unshaded.j3md");
//        		    gmat.getAdditionalRenderState().setWireframe(true);
//        		    gmat.getAdditionalRenderState().setLineWidth(1);
//        		    gmat.setColor("Color", ColorRGBA.Yellow);
//        		    gArrow.setMaterial(gmat);
//        		    rootNode.attachChild(gArrow);
//        		    gArrow.setLocalTranslation(r.get(i).pos);
        			
            		Material marker = new Material(assetManager, "Common/MatDefs/Misc/Unshaded.j3md");
                	marker.setColor("Color", new ColorRGBA(1f, (float)i / (float)r.size(), (float)(index) / rows.size(), 1f));
                	Geometry mg = new Geometry("VertexMarker", new Sphere(8, 8, 0.01f));
                	mg.setMaterial(marker);
                	mg.setLocalTranslation(r.get(i).pos);
                	rootNode.attachChild(mg);
                	
                	
            	}
        		System.out.println();
        		index++;
        	}
        	
//        	VertexBuffer vb = m.getBuffer(VertexBuffer.Type.Position);
//        	float[] vertices = BufferUtils.getFloatArray((FloatBuffer) vb.getData());
//        	
//        	Material marker = new Material(assetManager, "Common/MatDefs/Misc/Unshaded.j3md");
//        	marker.setColor("Color", ColorRGBA.Red);
//        	
//        	// For each vertex, create a small sphere and position it at the vertex's location
//            for (int i = 0; i < vertices.length; i += 3) {
//                Vector3f vertexPos = new Vector3f(vertices[i], vertices[i + 1], vertices[i + 2]);
//                Geometry mg = new Geometry("VertexMarker", new Sphere(8, 8, 0.01f)); // Small sphere
//                mg.setMaterial(marker);
//                mg.setLocalTranslation(vertexPos);
//                rootNode.attachChild(mg);
//            }
        }

        // Set up lighting
        AmbientLight ambient = new AmbientLight();
        ambient.setColor(ColorRGBA.White);
        rootNode.addLight(ambient);

        DirectionalLight sun = new DirectionalLight();
        sun.setDirection(new Vector3f(-1, -1, -1).normalizeLocal());
        sun.setColor(ColorRGBA.White);
        rootNode.addLight(sun);
        
        DirectionalLight backLight = new DirectionalLight();
        backLight.setDirection(new Vector3f(1, 1, 1).normalizeLocal());
        backLight.setColor(ColorRGBA.White);
        rootNode.addLight(backLight);

        cam.setLocation(new Vector3f(0, 0, 10));
        cam.lookAt(Vector3f.ZERO, Vector3f.UNIT_Y);
        flyCam.setMoveSpeed(10f);
        flyCam.setZoomSpeed(10f);
        flyCam.setDragToRotate(true);
        inputManager.setCursorVisible(true);
	}
	
	public Vector3f computeTangent(Vector3f gradF) {
		Vector3f gradient = gradF;
		
		Vector3f basis1 = gradient.normalize();
		
		Vector3f arbVec = new Vector3f(1f, 0f, 0f);
		if(Math.abs(basis1.dot(arbVec)) > 0.9d) {
			arbVec = new Vector3f(0f, 1f, 0f);
		}
		
		Vector3f basis2 = basis1.cross(arbVec).normalize();
		
		double coeff1 = gradient.dot(basis1);
		double coeff2 = gradient.dot(basis2);
		
		Vector3f tangent = basis2.mult((float) coeff1).subtract(basis1.mult((float) coeff2));
		
		return tangent;
	}
	
	public void resampleVertices(TreeMap<Double, List<Vertex>> treeMap) {
		double epsilon = GRID_WIDTH * 0.0001d;
		
		for(Entry<Double, List<Vertex>> rowEntry : treeMap.entrySet()) {
			
			Double key = rowEntry.getKey();
			List<Vertex> rowVertices = rowEntry.getValue();
			
			List<Vertex> selectedVertices = new ArrayList<>();
			selectedVertices.add(rowVertices.get(0));
			
			for(int i = 1; i < rowVertices.size(); i++) {
				double dist = rowVertices.get(i).pos.distance(selectedVertices.get(selectedVertices.size() - 1).pos);
				
				if(Math.abs(dist - GRID_WIDTH) <= epsilon) {
					selectedVertices.add(rowVertices.get(i));
				}else if(dist > GRID_WIDTH) {
					Vertex interVertex = interpolate(selectedVertices.get(selectedVertices.size() - 1), rowVertices.get(i), GRID_WIDTH);
					selectedVertices.add(interVertex);
					i--;
				}
			}
			
			treeMap.put(key, selectedVertices);
		}
	}
	
	public Vertex findOrCreateVertex(Vector3f pos, List<Vertex> verticesList) {
		for(Vertex v : verticesList) {
			if(v.pos.equals(pos)) {
				return v;
			}
		}
		
		Vertex newVertex = new Vertex();
		newVertex.pos = pos;
		verticesList.add(newVertex);
		return newVertex;
	}
	
	public void computeG(List<Vertex> verticesList, int maxIterations) {
		double learningRate = 0.01d;
		
		for(int iter = 0; iter < maxIterations; iter++) {
			
			for(Vertex v : verticesList) {
				v.gradF = computeGradientF(v);
				v.gradG = computeGradientG(v);
			}
			
			for(int i = 0; i < verticesList.size(); i++) {
				verticesList.get(i).rotatedGradF = computeTangent(verticesList.get(i).gradF);
				double dotProduct = verticesList.get(i).rotatedGradF.dot(verticesList.get(i).gradG);
				
				double error = dotProduct - 1;
				
				double gradientLocalObj = 2 * error * verticesList.get(i).rotatedGradF.length();
				
				verticesList.get(i).g -= learningRate * gradientLocalObj;
			}
		}
	}
	
	public Vector3f computeGradientG(Vertex v) {
		Vector3f gradient = new Vector3f();
		for(Vertex neighbor : v.neighbors) {
			Vector3f direction = neighbor.pos.subtract(v.pos);
			double dg = neighbor.g - v.g;
			Vector3f weightedDirection = direction.scaleAdd((float)dg, Vector3f.ZERO);
			gradient = gradient.add(weightedDirection);
		}
		
		return gradient.normalize();
	}
	
	public Vector3f computeGradientF(Vertex v) {
		Vector3f gradient = new Vector3f();
		for(Vertex neighbor : v.neighbors) {
			Vector3f direction = neighbor.pos.subtract(v.pos);
			double df = neighbor.f - v.f;
			Vector3f weightedDirection = direction.scaleAdd((float) df, Vector3f.ZERO);
			gradient = gradient.add(weightedDirection);
		}
		
		return gradient.normalize();
	}
	
//	public void updateTentativeDistance(Vertex current, Vertex neighbor) {
//		List<Vertex> knownNeighbors = new ArrayList<Vertex>();
//		knownNeighbors.add(current);
//		
//		knownNeighbors.addAll(neighbor.neighbors.stream().filter(v -> v.state == State.KNOWN).collect(Collectors.toList()));
//		
//		if (knownNeighbors.size() < 2) {
//	        return; // We need at least two known neighbors to update the distance.
//	    }
//
//	    // Sort the known neighbors by their distance.
//	    knownNeighbors.sort(Comparator.comparing(v -> v.f));
//
//	    Vertex a = knownNeighbors.get(0);
//	    Vertex b = knownNeighbors.get(1);
//
//	    double da = a.f;
//	    double db = b.f;
//
//	    double la = neighbor.pos.distance(a.pos);
//	    double lb = neighbor.pos.distance(b.pos);
//
//	    // Solve the quadratic equation derived from the Eikonal equation.
//	    double A = 2.0 / (la * la) + 2.0 / (lb * lb);
//	    double B = -2.0 * (da / (la * la) + db / (lb * lb));
//	    double C = (da * da) / (la * la) + (db * db) / (lb * lb) - 1.0;
//
//	    double discriminant = B * B - 4.0 * A * C;
//
//	    if (discriminant >= 0) {
//	        double T1 = (-B + Math.sqrt(discriminant)) / (2.0 * A);
//	        double T2 = (-B - Math.sqrt(discriminant)) / (2.0 * A);
//
//	        // Choose the smaller positive solution.
//	        double T = (T1 > 0 && T2 > 0) ? Math.min(T1, T2) : Math.max(T1, T2);
//
//	        if (T < neighbor.f) {
//	            neighbor.f = T;
//	        }
//	    }
//	}
//	
//	public void computeGeodesicDistanceFMM(List<Vertex> vertices, Vertex seed) {
//	    PriorityQueue<Vertex> queue = new PriorityQueue<>(Comparator.comparing(v -> v.f));
//
//	    for (Vertex v : vertices) {
//	        v.f = Double.POSITIVE_INFINITY;
//	        v.state = State.FAR;
//	    }
//
//	    seed.f = 0;
//	    seed.state = State.KNOWN;
//	    
//	    System.out.println("Total vertices: " + vertices.size());
//	    System.out.println("Seed vertex: " + seed.pos + " with init dist: " + seed.f);
//	    System.out.println("Number of seed neighbors: " + seed.neighbors.size());
//
//	    for (Vertex neighbor : seed.neighbors) {
//	        neighbor.f = seed.pos.distance(neighbor.pos);
//	        neighbor.state = State.TRIAL;
//	        queue.add(neighbor);
////	        System.out.println("Seed neighbors: " + neighbor.pos + " with dist: " + neighbor.f + " and state: " + neighbor.state);
//	    }
//
//	    while (!queue.isEmpty()) {
//	        Vertex current = queue.poll();
//	        current.state = State.KNOWN;
//	        long kc = current.neighbors.stream().filter(v -> v.state == State.KNOWN).count();
////	        System.out.println("Processing vertex: " + current.pos + " with dist: " + current.f + " and state: " + current.state + " knowns: " + kc);
//	        
//	        for (Vertex neighbor : current.neighbors) {
//	        	double oldDist = neighbor.f;
//	            if (neighbor.state == State.FAR) {
//	                updateTentativeDistance(current, neighbor);
//	                neighbor.state = State.TRIAL;
//	                queue.add(neighbor);
////	                System.out.println("Neighbor: " + neighbor.pos + " with old dist: " + oldDist + " and new dist: " + neighbor.f + " and state: " + neighbor.state);
//	            }else {
////	            	System.out.println("Neighbor (not updated): " + neighbor.pos + " with dist: " + neighbor.f + " and state: " + neighbor.state);
//	            }
//	        }
////	        System.out.println();
//	        
////	        for(Vertex vt : queue.stream().toList()) {
////	        	System.out.print(vt.pos.toString() + " " + vt.f + ") ");
////	        }
////	        System.out.println();
//	    }
//	    
//	    long knownCount = vertices.stream().filter(v -> v.state == State.KNOWN).count();
//	    long trialCount = vertices.stream().filter(v -> v.state == State.TRIAL).count();
//	    long farCount = vertices.stream().filter(v -> v.state == State.FAR).count();
////	    System.out.println("Vertices with state KNOWN: " + knownCount);
////	    System.out.println("Vertices with state TRIAL: " + trialCount);
////	    System.out.println("Vertices with state FAR: " + farCount);
//
//	}
	
	public void computeGeodesicDistance(List<Vertex> vertices, Vertex seed) {
		PriorityQueue<Vertex> queue = new PriorityQueue<>(Comparator.comparing(v -> v.f));
		
		seed.f = 0;
		queue.add(seed);
		
		while(!queue.isEmpty()) {
			Vertex current = queue.poll();
			current.finalized = true;
			
			for(Vertex neighbor : current.neighbors) {
				if(!neighbor.finalized) {
					double newDist = current.f + current.pos.distance(neighbor.pos);
					
					if(newDist < neighbor.f) {
						neighbor.f = newDist;
						queue.remove(neighbor);
						queue.add(neighbor);
					}
				}
			}
		}
	}
	
	public static boolean areClose(double a, double b) {
		return Math.abs(a - b) < THRESHOLD;
	}
	
	public static TreeMap<Double, List<Vertex>> groupByRows(List<Vertex> vertices) {
		TreeMap<Double, List<Vertex>> rows = new TreeMap<>();
		
		for(Vertex vertex : vertices) {
			boolean added = false;
			
			for(Double key : rows.keySet()) {
				if(areClose(vertex.f, key)) {
					rows.get(key).add(vertex);
					added = true;
					break;
				}
			}
			
			if(!added) {
				List<Vertex> newRow = new ArrayList<>();
				newRow.add(vertex);
				rows.put(vertex.f, newRow);
			}
		}
		
		return rows;
	}
	
	public static Vertex findClosestVertex(Vertex input, List<Vertex> vertices) {
		
		if(vertices == null || vertices.isEmpty()) {
			return null;
		}
		
		Vertex closest = vertices.get(0);
		double minDist = input.pos.distance(closest.pos);
		for(Vertex vertex : vertices) {
			
			double curDist = input.pos.distance(vertex.pos);
			if(curDist < minDist) {
				minDist = curDist;
				closest = vertex;
			}
		}
		
		return closest;
	}
	
	public static Vertex interpolate(Vertex start, Vertex end, double targetDist) {
		Vector3f direction = end.pos.subtract(start.pos);
		
		direction.normalizeLocal();
		
		direction.multLocal((float) targetDist);
		
		Vector3f interpolatedPos = start.pos.add(direction);
		
		double totalDist = start.pos.distance(end.pos);
		double weight = targetDist / totalDist;
		
		double interpolatedF = start.f + weight * (end.f - start.f);
		double interpolatedG = start.g + weight * (end.g - start.g);
		
		Vertex interpolatedVertex = new Vertex();
		interpolatedVertex.pos = interpolatedPos;
		interpolatedVertex.f = interpolatedF;
		interpolatedVertex.g = interpolatedG;
		
		return interpolatedVertex;
	}
}