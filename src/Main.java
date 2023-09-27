import java.awt.DisplayMode;
import java.awt.GraphicsDevice;
import java.awt.GraphicsEnvironment;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.TreeMap;

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
	
	public static final double THRESHOLD = 1e-3;
	public static final double GRID_WIDTH = 0.5d;
	
	private List<Vertex> verticesList;

	public static void main(String[] args) {
		Main main = new Main();
		
		GraphicsDevice device = GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice();
		DisplayMode[] modes = device.getDisplayModes();
		int i = 0;
		AppSettings settings = new AppSettings(true);
		settings.setResolution(1080, 720);
		settings.setFrequency(modes[i].getRefreshRate());
		settings.setBitsPerPixel(modes[i].getBitDepth());
		settings.setVSync(true);
		settings.setGammaCorrection(true);
		main.setSettings(settings);
		main.setShowSettings(false);
		
		main.setDisplayFps(true);
		
		main.start();
	}
	
	@Override
	public void simpleInitApp() {
		Spatial mesh = assetManager.loadModel("sphere2.obj");
		Material mat = new Material(assetManager, "Common/MatDefs/Light/Lighting.j3md");
        mat.setColor("Diffuse", ColorRGBA.Gray);
        mat.setColor("Ambient", ColorRGBA.White.mult(0.1f));
        mat.setBoolean("UseMaterialColors", true);
        mesh.setMaterial(mat);
        rootNode.attachChild(mesh);
        mesh.setLocalTranslation(-3f, 0f, 0f);
        
        verticesList = new ArrayList<Vertex>();
        
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
        	
        	Vertex seed = findClosestVertex(new Vector3f(0f, -1f, 0f), verticesList);
        	
        	computeGeodesicDistance(verticesList, seed);
        	
        	TreeMap<Double, List<Vertex>> rows = groupByRows(verticesList);
        	
        	List<Vertex> boundary = findLongestGeodesic(verticesList, seed);
        	
        	for(Vertex v : verticesList) {
        		for(Vertex b : boundary) {
        			if(v.equals(b)) {
        				v.isBoundary = true;
        				v.g = 0d;
        			}
        		}
        	}
        	
        	computeG(verticesList, 200);
        	
//        	resampleVertices(rows);
        	
        	int index = 0;
        	for(Map.Entry<Double, List<Vertex>> entry : rows.entrySet()) {
        		List<Vertex> r = entry.getValue();
        		
        		if(index < 11 || index > 11) {
        			index++;
        			continue;
        		}
        		
        		System.out.println("row size: " + r.size());
        		for(int i = 0; i < r.size(); i++) {
        			System.out.println("g value: " + r.get(i).g + " f value: " + r.get(i).f + " grad g: " + r.get(i).gradG);
//        			Vector3f arrowVec = new Vector3f(r.get(i).rotatedGradF);
//        			Arrow arrow = new Arrow(arrowVec);
//        			Geometry gArrow = new Geometry("f vector", arrow);
//        			Material gmat = new Material(assetManager, "Common/MatDefs/Misc/Unshaded.j3md");
//        		    gmat.getAdditionalRenderState().setWireframe(true);
//        		    gmat.getAdditionalRenderState().setLineWidth(4);
//        		    gmat.setColor("Color", ColorRGBA.Green);
//        		    gArrow.setMaterial(gmat);
//        		    rootNode.attachChild(gArrow);
//        		    gArrow.setLocalTranslation(r.get(i).pos);
//        		    
        		    Vector3f arrowVec2 = new Vector3f(r.get(i).gradG);
        		    Arrow arrow2 = new Arrow(arrowVec2);
        		    Geometry gArrow2 = new Geometry("grad f vector", arrow2);
        		    Material gmat2 = new Material(assetManager, "Common/MatDefs/Misc/Unshaded.j3md");
        		    gmat2.getAdditionalRenderState().setWireframe(true);
        		    gmat2.getAdditionalRenderState().setLineWidth(4);
        		    gmat2.setColor("Color", ColorRGBA.Red);
        		    gArrow2.setMaterial(gmat2);
        		    rootNode.attachChild(gArrow2);
        		    gArrow2.setLocalTranslation(r.get(i).pos);
        		    
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
        			
        		    float hueForRow = (float) (r.get(i).f);
        		    float hueForVertex = ((float)(r.get(i).g));
        		    float hue = (hueForRow + hueForVertex) / 2f;
        		    ColorRGBA color = hslToRgb(hue, 0.5f, 0.5f);
        		    
            		Material marker = new Material(assetManager, "Common/MatDefs/Misc/Unshaded.j3md");
                	marker.setColor("Color", color);
                	Geometry mg = new Geometry("VertexMarker", new Sphere(8, 8, 0.01f));
                	mg.setMaterial(marker);
                	mg.setLocalTranslation(r.get(i).pos);
                	rootNode.attachChild(mg);
                	
                	Material marker2 = new Material(assetManager, "Common/MatDefs/Misc/Unshaded.j3md");
                	marker2 = new Material(assetManager, "Common/MatDefs/Misc/Unshaded.j3md");
                	marker2.setColor("Color", color);
                	Geometry mg2 = new Geometry("VertexMarker", new Sphere(8, 8, 0.01f));
                	mg2.setMaterial(marker2);
//                	Vector3f vec = new Vector3f(r.get(i).pos.x, r.get(i).pos.y, r.get(i).pos.z);
                	Vector3f vec = new Vector3f((float)(r.get(i).g), (float)r.get(i).f, 0f);
                	vec.addLocal(3f, 0f, 0f);
                	mg2.setLocalTranslation(vec);
                	rootNode.attachChild(mg2);
            	}
        		System.out.println();
        		index++;
        	}
        	
        	System.out.println(boundary.size());
        	for(Vertex v : boundary) {
        		Material marker = new Material(assetManager, "Common/MatDefs/Misc/Unshaded.j3md");
            	marker.setColor("Color", new ColorRGBA(0f, 0f, 1f, 1f));
            	Geometry mg = new Geometry("VertexMarker", new Sphere(8, 8, 0.01f));
            	mg.setMaterial(marker);
            	mg.setLocalTranslation(v.pos);
            	rootNode.attachChild(mg);
        	}
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
	
	public Vector3f computeTangent(Vertex v) {
		// Step 1: Normalize gradF
	    Vector3f normalizedGradF = v.gradF.normalize();

	    // Step 2: Compute the normal at the vertex (since it's a sphere)
	    Vector3f normal = v.pos.normalize();

	    // Step 3: Compute a vector orthogonal to gradF in the tangent plane
	    Vector3f tangent = normalizedGradF.cross(normal);

	    // Step 4: Normalize the tangent vector
	    tangent.normalizeLocal();

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
		int n = verticesList.size();
		double[][] H = new double[n][n];
		for(int i = 0; i < n; i++) {
			H[i][i] = 1.0d;
		}
		
		for(int iter = 0; iter < maxIterations; iter++) {
			double[] localGradients = new double[n];
			double[] p = new double[n]; // search direction
			
			double sumOfErrors = 0d;
			
			for(int i = 0; i < n; i++) {
				Vertex v = verticesList.get(i);
				v.gradG = computeGradientG(v);
				v.gradF = computeGradientF(v);
				v.rotatedGradF = computeTangent(v);
				
				double error = v.rotatedGradF.dot(v.gradG) - 1d;
				sumOfErrors += error * error;
				
				localGradients[i] = 2d * error;
			}
			
			double meanError = sumOfErrors / (double)n;
			
			System.out.println("epoch: " + iter + " sum of errors: " + sumOfErrors + " mean error: " + meanError);
			
			for(int i = 0; i < n; i++) {
				for(int j = 0; j < n; j++) {
					p[i] -= H[i][j]	* localGradients[j];
				}
			}
			
			
			
			double[] s = new double[n];
			double[] y = new double[n];
			double[] newLocalGradients = new double[n];
			for(int i = 0; i < n; i++) {
				Vertex v = verticesList.get(i);
				if(!v.isBoundary) {
					double alpha = lineSearch(i, p, localGradients, 1.0d, 0.01d, 0.5d); // should use line search to find alpha
					s[i] = alpha * p[i];
					v.g += s[i];
					v.gradG = computeGradientG(v);
					v.gradF = computeGradientF(v);
					v.rotatedGradF = computeTangent(v);
					newLocalGradients[i] = 2d * (v.rotatedGradF.dot(v.gradG) - 1d);
					y[i] = newLocalGradients[i] - localGradients[i];
				}
			}
			
			double yTs = 0;
			for(int i = 0; i < n; i++) {
				yTs += y[i] * s[i];
			}
			
			double[][] term1 = new double[n][n];
			double[][] term2 = new double[n][n];
			double[][] term3 = new double[n][n];
			for(int i = 0; i < n; i++) {
				for(int j = 0; j < n; j++) {
					term1[i][j] = s[i] * s[j] / yTs;
					term2[i][j] = s[i] * y[j] / yTs;
					term3[i][j] = y[i] * s[j] / yTs;
				}
			}
			
			double[][] H_new = new double[n][n];
			for(int i = 0; i < n; i++) {
				for(int j = 0; j < n; j++) {
					H_new[i][j] = H[i][j] - H[i][j] * term2[i][j] - term3[i][j] * H[i][j] + term1[i][j];
				}
			}
			
			H = H_new;
			
			localGradients = newLocalGradients;
		}
	}
	
	public double lineSearch(int vertexIndex, double[] p, double[] gradient, double alphaStart, double c, double rho) {
		Vertex v = verticesList.get(vertexIndex);
		double alpha = alphaStart;
		double currentObjective = computeObjective(v, v.g);
		while(computeObjective(v, v.g + alpha * p[vertexIndex]) > currentObjective + c * alpha * gradient[vertexIndex] * p[vertexIndex]) {
			alpha *= rho;
		}
		
		return alpha;
	}
	
	private double computeObjective(Vertex v, double gValue) {
		double tempG = v.g;
		v.g = gValue;
		v.gradG = computeGradientG(v);
		v.gradF = computeGradientF(v);
		v.rotatedGradF = computeTangent(v);
		double error = v.rotatedGradF.dot(v.gradG) - 1d;
		v.g = tempG;
		return error * error;
	}
	
//	public void computeG(List<Vertex> verticesList, int maxIterations) {
//		Vector3f[] gradients = new Vector3f[verticesList.size()];
//		Vector3f[] directions = new Vector3f[verticesList.size()];
//		double[] g_values = new double[verticesList.size()];
//		
//		for(int i = 0; i < verticesList.size(); i++) {
//			Vertex v = verticesList.get(i);
//			gradients[i] = computeGradientG(v);
//			directions[i] = gradients[i].negate();
//			g_values[i] = v.g;
//		}
//		
//		for(int iter = 0; iter < maxIterations; iter++) {
//			
//			double alpha = 0.01d;
//			
//			for(int i = 0; i < verticesList.size(); i++) {
//				Vertex v = verticesList.get(i);
//				if(v.isBoundary) {
//					continue;
//				}
//				
//				g_values[i] += alpha * directions[i].length();
//				gradients[i] = computeGradientG(v);
//			}
//			
//			for(int i = 0; i < verticesList.size(); i++) {
//				Vertex v = verticesList.get(i);
//				if(v.isBoundary) {
//					continue;
//				}
//				
//				double beta = 0.1d;
//				
//				directions[i] = gradients[i].negate().add(directions[i].mult((float) beta));
//			}
//			
//			for(int i = 0; i < verticesList.size(); i++) {
//				Vertex v = verticesList.get(i);
//				v.g = g_values[i];
//			}
//		}
//		
//		for(int i = 0; i < verticesList.size(); i++) {
//			verticesList.get(i).gradG = gradients[i];
//		}
//	}
	
//	public void computeG(List<Vertex> verticesList, int maxIterations) {
//		double learningRate = 0.01d;
//		double momentum = 0.9d;
//		double prevAverageError = Double.MAX_VALUE;
//		
//		for(int iter = 0; iter < maxIterations; iter++) {
//			double sums = 0d;
//			double averageError = 0d;
//			
//			for(Vertex v : verticesList) {
//				
//				if(v.isBoundary) {
//					v.g = 0d;
//					continue;
//				}
//				
//				v.gradF = computeGradientF(v);
//				v.gradG = computeGradientG(v);
//				v.rotatedGradF = computeTangent(v);
//				
//				double dotProduct = v.rotatedGradF.dot(v.gradG);
//				double error = dotProduct - 1d;
//				double squaredError = error * error;
//				
//				// compute local gradient of the objective function
//				Vector3f gradientLocalObjective = v.rotatedGradF.mult((float) (2d * error));
//				
//				double gradJProj = 2d * error * gradientLocalObjective.dot(v.rotatedGradF);
//				
//				v.velocity = momentum * v.velocity -learningRate * gradJProj;
//				v.g += v.velocity;
//				
//				sums += squaredError;
//			}
//			
//			averageError = sums / (double)verticesList.size();
//			
//			System.out.println("epoch: " + iter + " sum: " + sums + " average error: " + averageError);
//			
////			if(Math.abs(prevAverageError - averageError) < 1e-10) {
////				System.out.println("converged at iteration " + iter);
////				break;
////			}
//			prevAverageError = averageError;
//		}
//	}
	
//	public void computeG(List<Vertex> verticesList, int maxIterations) {
//		double learningRate = 0.001d;
//		
//		for(int iter = 0; iter < maxIterations; iter++) {
//			double sums = 0d;
//			double averageError = 0d;
//			
//			for(Vertex v : verticesList) {
//				
//				if(v.isBoundary) {
//					v.g = 0d;
//					continue;
//				}
//				
//				v.gradF = computeGradientF(v);
//				v.gradG = computeGradientG(v);
//				v.rotatedGradF = computeTangent(v);
//				
//				double dotProduct = v.rotatedGradF.dot(v.gradG);
//				
//				double error = dotProduct - 1d;
//				
//				double squaredError = error * error;
//				
//				Vector3f gradientLocalObjective = v.rotatedGradF.mult((float) (2d * error));
//				
//				v.gradG = v.gradG.subtract(gradientLocalObjective.mult((float)learningRate));
//				
//				double deltaG = -learningRate * v.gradG.length();
//				v.g += deltaG;
//				
//				sums += squaredError;
//			}
//			
//			averageError = sums / (double)verticesList.size();
//			
//			System.out.println("epoch: " + iter + " sum: " + sums + " average error: " + averageError);
//		}
//	}
	
//	public void computeG(List<Vertex> verticesList, int maxIterations) {
//		double learningRate = 0.01d;
//		double momentum = 0.9d;
//		
//		for(int iter = 0; iter < maxIterations; iter++) {
////			System.out.println("Initial g values: " + verticesList.stream().map(v -> v.g).collect(Collectors.toList()));
//			
//			double sums = 0d;
//			double averageError = 0d;
//			double[] localGradients = new double[verticesList.size()];
//			
//			for(int i = 0; i < verticesList.size(); i++) {
//				Vertex v = verticesList.get(i);
//				v.gradF = computeGradientF(v);
//				v.gradG = computeGradientG(v);
//				v.rotatedGradF = computeTangent(v);
//				
//				double dotProduct = v.rotatedGradF.dot(v.gradG);
//				
//				double error = dotProduct - 1d;
//				
//				double squaredError = error * error;
//				
//				double gradientLocalObjective = 2d * error;
//				
//				localGradients[i] = gradientLocalObjective;
//				
//				sums += squaredError;
////				System.out.println("v.g: " + v.g + " learningrate: " + learningRate + " glo: " + gradientLocalObjective);
//			}
//			
//			for(int i = 0; i < verticesList.size(); i++) {
//				Vertex v = verticesList.get(i);
//				if(!v.isBoundary) {
//					v.velocity = momentum * v.velocity + learningRate * localGradients[i];
//					v.g -= v.velocity;
//				}
//			}
//			
//			averageError = sums / (double)verticesList.size();
//			
////			System.out.println("Updated g values: " + verticesList.stream().map(v -> v.g).collect(Collectors.toList()));
//			System.out.println("epoch: " + iter + " sum: " + sums + " average error: " + averageError);
//		}
//	}
	
//	public void computeG(List<Vertex> verticesList, int maxIterations) {
//		double learningRate = 0.01d;
//		double momentum = 0.9d;
//		
//		for(int iter = 0; iter < maxIterations; iter++) {
////			System.out.println("Initial g values: " + verticesList.stream().map(v -> v.g).collect(Collectors.toList()));
//			
//			double sums = 0d;
//			double averageError = 0d;
//			
//			for(Vertex v : verticesList) {
//				
//				if(v.isBoundary) {
//					v.g = 0d;
//					continue;
//				}
//				
//				v.gradF = computeGradientF(v);
//				v.gradG = computeGradientG(v);
//				v.rotatedGradF = computeTangent(v);
//				
//				double dotProduct = v.rotatedGradF.dot(v.gradG);
//				
//				double error = dotProduct - 1d;
////				System.out.println("dot: " + dotProduct + " " + " error: " + error);
//				
//				double squaredError = error * error;
//				
//				double gradientLocalObjective = 2d * error * v.rotatedGradF.length();
//				
//				sums += squaredError;
//				
//				v.velocity = momentum * v.velocity + learningRate * gradientLocalObjective;
//				v.g -= v.velocity;
////				System.out.println("v.g: " + v.g + " learningrate: " + learningRate + " glo: " + gradientLocalObjective);
//				
//			}
//			
//			averageError = sums / (double)verticesList.size();
//			
////			System.out.println("Updated g values: " + verticesList.stream().map(v -> v.g).collect(Collectors.toList()));
//			System.out.println("epoch: " + iter + " sum: " + sums + " average error: " + averageError);
//		}
//	}
	
	public Vector3f computeGradientG(Vertex v) {
		Vector3f gradient = new Vector3f();
		for(Vertex neighbor : v.neighbors) {
			Vector3f direction = neighbor.pos.subtract(v.pos);
			double dg = neighbor.g - v.g;
			Vector3f weightedDirection = direction.scaleAdd((float)dg, Vector3f.ZERO);
			gradient = gradient.add(weightedDirection);
		}
		
		return gradient.normalizeLocal();
	}
	
	public Vector3f computeGradientF(Vertex v) {
		Vector3f gradient = new Vector3f();
		for(Vertex neighbor : v.neighbors) {
			Vector3f direction = neighbor.pos.subtract(v.pos);
			double df = neighbor.f - v.f;
			Vector3f weightedDirection = direction.scaleAdd((float) df, Vector3f.ZERO);
			gradient = gradient.add(weightedDirection);
		}
		
		return gradient.normalizeLocal();
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
	
	public List<Vertex> findLongestGeodesic(List<Vertex> vertices, Vertex seed) {
	    // Step 1: Compute shortest paths using Dijkstra's algorithm
	    PriorityQueue<Vertex> queue = new PriorityQueue<>(Comparator.comparingDouble(v -> v.distance));
	    seed.distance = 0;
	    queue.add(seed);

	    while (!queue.isEmpty()) {
	        Vertex current = queue.poll();

	        if (current.distance == Double.MAX_VALUE) {
	            break;  // All remaining vertices are unreachable from the seed
	        }

	        for (Vertex neighbor : current.neighbors) {
	            double newDist = current.distance + current.pos.distance(neighbor.pos);
	            if (newDist < neighbor.distance) {
	                queue.remove(neighbor);  // Remove the old instance from the queue
	                neighbor.distance = newDist;
	                neighbor.predecessor = current;
	                queue.add(neighbor);  // Add the updated instance to the queue
	            }
	        }
	    }

	    // Step 2: Identify the vertex with the maximum geodesic distance
	    Vertex farthest = seed;
	    for (Vertex v : vertices) {
	        if (v.distance > farthest.distance) {
	            farthest = v;
	        }
	    }

	    // Step 3: Reconstruct the path from the farthest vertex to the seed
	    List<Vertex> path = new ArrayList<>();
	    while (farthest != null) {
	        path.add(farthest);
	        farthest = farthest.predecessor;
	    }
	    Collections.reverse(path);  // so the path starts from the seed

	    return path;
	}

 	
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
	
	public static Vertex findClosestVertex(Vector3f input, List<Vertex> vertices) {
		
		if(vertices == null || vertices.isEmpty()) {
			return null;
		}
		
		Vertex closest = vertices.get(0);
		double minDist = input.distance(closest.pos);
		for(Vertex vertex : vertices) {
			
			double curDist = input.distance(vertex.pos);
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
	
	public ColorRGBA hslToRgb(float h, float s, float l){
	    float r, g, b;

	    if(s == 0){
	        r = g = b = l; // achromatic
	    }else{
	        float q = l < 0.5 ? l * (1 + s) : l + s - l * s;
	        float p = 2 * l - q;
	        r = hueToRgb(p, q, h + 1/3f);
	        g = hueToRgb(p, q, h);
	        b = hueToRgb(p, q, h - 1/3f);
	    }
	    return new ColorRGBA(r, g, b, 1f);
	}

	private float hueToRgb(float p, float q, float t){
	    if(t < 0) t += 1;
	    if(t > 1) t -= 1;
	    if(t < 1/6f) return p + (q - p) * 6 * t;
	    if(t < 1/2f) return q;
	    if(t < 2/3f) return p + (q - p) * (2/3f - t) * 6;
	    return p;
	}
}