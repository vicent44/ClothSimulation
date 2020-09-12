﻿using System.Collections;
using System.Collections.Generic;
using UnityEngine;


[RequireComponent(typeof(MeshFilter))]
public class MeshGenerator : MonoBehaviour
{
    //Create the mesh in unity
    Mesh mesh; 
    Vector3[] vertices;
    int[] triangles;
    //Parameter to set the size of the grid
    public int gridSizeNew;
    //Variable to choose if we want to draw the lines between particles
    [SerializeField] bool drawSprings = true;
    //Setting the cloth density
    [SerializeField] float clothDensity = 150f;
    //Wind direction, in normal vector
    [SerializeField] Vector3 windDirection=new Vector3(0,0,1);
    //Wind module
    [SerializeField] float windModule = 10f;

    //Setting all the constants for the springs
    [SerializeField] float elasticConstantStructural = 20f;
    [SerializeField] float elasticConstantShear = 20f;
    [SerializeField] float elasticConstantBend = 20f;
    [SerializeField] float dampingConstant = 0.7f;
    //Setting all the constants of friction and dissipation for cloth and plane
    [SerializeField] float frictionConstPlane = 0.7f;
    [SerializeField] float dissipationConstPlane = 0.1f;
    [SerializeField] float frictionConstCloth = 0.7f;
    [SerializeField] float dissipationConstCloth = 0.1f;

    Simulate simulator;
    Hashing hash;

    float timePassed = 0.0f;
    public float deltaTimeStep = 0.02f;

    //List with all the objects: particles, springs, triangles.
    List<Particles> _particles;
    List<Springs> _springs;
    List<Triangles> _triangles;

    Vector3 winddirectiondensity;
    //List of gameobjects to be able to select some vertex in the mesh (gameobjects)
    [SerializeField] GameObject sphereControl;
    List<GameObject> _sphere;
    [SerializeField] float sphereScale = 0.1f;
    //Set the transform plane to know where to do the plane collision
    [SerializeField] Transform plane;
    [SerializeField] Transform secondPlane;
    [SerializeField] Vector3 normalSecondPlane = new Vector3(0f, 0f, 1f);

    //Mouse Drag Variables
    Vector3 screenPoint;
    Vector3 offset;
    Ray ray;
    RaycastHit hit;

    public bool isPaused = true;

    void Start()
    {
        //Setting the size of the mesh, specially the density
        //int gridSize = gridSizeNew + (gridSizeNew - 1);
        int gridSize = gridSizeNew + gridSizeNew * (gridSizeNew - 1);
        
        //Initialice the list of particles, springs and triangles.
        _particles = new List<Particles>();
        _springs = new List<Springs>();
        _triangles = new List<Triangles>();
        _sphere = new List<GameObject>();

        //To be able to change the size of the sphere
        sphereControl.transform.localScale = new Vector3(sphereScale, sphereScale, sphereScale);

        //In the editor to have a parent of all the particles
        var particleParent = new GameObject("particleParent");

        //Calcule the area of the cloth
        float area = 0.2f * (gridSize - 1) * 0.2f * (gridSize - 1);

        //Calcule the total mass of the cloth with the density
        float massTotal = clothDensity * area;
        //Calcule the mass of each particle
        float mass = massTotal / (gridSize * gridSize);

        //Instantiate the position of the particles and the gameobject.
        float jP = 0;
        float iP = 0;
        int pi = 0;
        int pj = 0;
        for(int i = 0; i < gridSize; i++)
        {
            for(int j = 0; j < gridSize; j++)
            {
                jP = j * 0.2f;
                iP = i * 0.2f;
                
                var newP = ParticlesBehaviour.Create(new Vector3(jP,0.0f,iP), mass, pi, pj, sphereControl);
                newP.transform.SetParent(particleParent.transform, false);
                newP.transform.localPosition = new Vector3(jP,0.0f,iP);
                newP.particles.Position = newP.transform.position;
                _particles.Add(newP.particles);
                pi++;
            }
            pj++;
        }
        
        //Springs
        //Structural Horizontal - Right
        for(int j = 0; j < gridSize; j++)
        {
            for(int i = 0; i < gridSize - 1; i++)
            {
                var index = j * gridSize + i;
                var c = _particles[index];
                var right = _particles[index + 1];
                var re = new Springs(c, right, elasticConstantStructural, dampingConstant, 1);
                _springs.Add(re);
            }
        }

        //Structural Vertical - Bot
        for(int j = 0; j < gridSize - 1; j++)
        {
            for(int i = 0; i < gridSize; i++)
            {
                var index_1 = j * gridSize + i;
                var index_2 = (j + 1) * gridSize + i;
                var c = _particles[index_1];
                var right = _particles[index_2];
                var re = new Springs(c, right, elasticConstantStructural, dampingConstant, 1);
                _springs.Add(re);
            }
        }

        //Shear - \ - Down
        for(int j = 0; j < gridSize - 1; j++)
        {
            for(int i = 0; i < gridSize - 1; i++)
            {
                var index_1 = j * gridSize + i;
                var index_2 = (j + 1) * gridSize + i + 1;
                var c = _particles[index_1];
                var bot = _particles[index_2];
                var re = new Springs(c, bot, elasticConstantShear, dampingConstant, 2);
                _springs.Add(re);
            }
        }

        //Shear - / - Down
        for(int j = 0; j < gridSize - 1; j++)
        {
            for(int i = 1; i < gridSize; i++)
            {
                var index_1 = j * gridSize + i;
                var index_2 = (j + 1) * gridSize + i - 1;
                var c = _particles[index_1];
                var top = _particles[index_2];
                var re = new Springs(c, top, elasticConstantShear, dampingConstant, 2);
                _springs.Add(re);
            }
        }

        //Bend - Left
        for(int j = 0; j < gridSize; j++)
        {
            for(int i = 0; i < gridSize - 2; i++)
            {
                var index_1 = j * gridSize + i;
                var index_2 = j * gridSize + i + 2;
                var c = _particles[index_1];
                var left = _particles[index_2];
                var re = new Springs(c, left, elasticConstantBend, dampingConstant, 3);
                _springs.Add(re);
            }
        }

        //Bend - Down
        for(int j = 0; j < gridSize - 2; j++)
        {
            for(int i = 0; i < gridSize; i++)
            {
                var index_1 = j * gridSize + i;
                var index_2 = (j + 2) * gridSize + i;
                var c = _particles[index_1];
                var left = _particles[index_2];
                var re = new Springs(c, left, elasticConstantBend, dampingConstant, 3);
                _springs.Add(re);
            }
        }
        
        //Fix two particles
        _particles[0].isActive = false;
        _particles[gridSize - 1].isActive = false;

        Vector3.Normalize(normalSecondPlane);

        //Calcule a part of the wind force and call the simulation script to incitialize it with all the needed information
        winddirectiondensity = windDirection * windModule * clothDensity;
        simulator = new Simulate(_particles, _springs, _triangles, winddirectiondensity, plane, secondPlane, normalSecondPlane, gridSize, frictionConstPlane, dissipationConstPlane, frictionConstCloth, dissipationConstCloth, drawSprings);

        //Creating the mess
        mesh = new Mesh();
        GetComponent<MeshFilter>().mesh = mesh;
        CreateShape();
        UpdateMesh();
    }

    void FixedUpdate()
    {
        
        /*if(Input.GetKey("left"))
        {
            Debug.Log("hola");
            _particles[0].isActive = true;
            _particles[0].AddPosition(new Vector3(0.01f, 0f, 0f));
            //_particles[0].ResetResultantForce();
            //_particles[0].Velocity = Vector3.zero;
            
            //_particles[0].ResetResultantForce();
        }
        if(Input.GetKey("right"))
        {
            _particles[0].isActive = true;
            _particles[0].Position -= new Vector3(0.01f, 0f, 0f);
        }
        if(Input.GetKey("up"))
        {
            _particles[0].Position += new Vector3(0f, 0f, 0.01f);
        }
        if(Input.GetKey("down"))
        {
            _particles[0].Position -= new Vector3(0.1f, 0f, 0.01f);
        }
        */
        //_particles[0].isActive = false;
        /*timePassed += Time.fixedDeltaTime;
        if(timePassed >= deltaTimeStep) timePassed = 0.0f;
        if(!isPaused && timePassed == 0.0f)
        {
            int gridSize = gridSizeNew + gridSizeNew * (gridSizeNew - 1);
            simulator.Update(deltaTimeStep, _triangles);
            UpdateMesh();
            Debug.Log(deltaTimeStep);
        }*/
        /*if(Input.GetKey("up"))
        {
            Debug.Log("hola");
        }*/
        //Collision constrain

        /*hash = new Hashing(gridSize, 1.0f/(float)gridSize, 1723);
        SPHash[] spHash = new SPHash[1723];

        for(int v = 0; v < _particles.Count; v++)
        {
            int has = Mathf.Abs(hash.Hash(_particles[v].Position));

            if(spHash[has].indices == null) spHash[has].indices = new List<int>();
            spHash[has].indices.Add(_particles[v].I);
        }

        for(int t = 0; t < _triangles.Count; t++)
        {
            var tri = _triangles[t];
            var p0 = _particles[tri.indexTriA].Position;
            var p1 = _particles[tri.indexTriB].Position;
            var p2 = _particles[tri.indexTriC].Position;

            float w0 = 1f/_particles[tri.indexTriA].Mass;
            float w1 = 1f/_particles[tri.indexTriB].Mass;
            float w2 = 1f/_particles[tri.indexTriC].Mass;

            List<int> hashes = hash.TriangleBoundingBoxHashes(p0, p1, p2);

            for(int h = 0; h < hashes.Count; h++)
            {
                if(spHash[h].indices != null)
                {
                    for(int sph = 0; sph < spHash[h].indices.Count; sph++)
                    {
                        int idx = spHash[h].indices[sph];
                        if(idx != _particles[tri.indexTriA].I && idx != _particles[tri.indexTriB].I && idx != _particles[tri.indexTriC].I)
                        {
                            Vector3 p = _particles[idx].Position;
                            float w = _particles[idx].Mass;

                            Vector3 corr, corr0, corr1, corr2;
                            if(intersec.TrianglePointDistance(
                                p, w,
                                p0, w0,
                                p1, w1,
                                p2, w2,
                                0.02f, 1f, 0.0f,
                                out corr, out corr0, out corr1, out corr2))
                            {
                                Debug.Log("Collition");
                                //_particles[idx].Position += corr;
                                //_particles[tri.indexTriA].Position += corr0;
                                //_particles[tri.indexTriB].Position += corr1;
                                //_particles[tri.indexTriC].Position += corr2;
                            }
                        }
                    }
                }
            }
        }*/

        //_particles[0].isActive = false;
        /*for(int i = 1; i < gridSize; i++)
        {
            for(int j = 1; j <= 2 * (gridSize - 1); j++)
            {
                pos1 = i * (2 * (gridSize - 1)) + j;
                pos2 = (i + 1) * (2 * (gridSize - 1)) + j;
                if(_triangles[pos1].numTri < _triangles[pos1 + 3].numTri)
            }
        }
        */


        //FindMaxAge(_springs);

        /*if(Input.GetMouseButtonDown(0))
        //if(Event.current.type == EventType.MouseDrag)
        {
            ray = Camera.main.ScreenPointToRay(Input.mousePosition);
            if(Physics.Raycast(ray, out hit))
            {
                Debug.Log(hit.point);
            }
        }*/


        /*for(int i = 0; i < _particles.Count; i++)
        {
            _esfere[i].transform.position = _particles[i].Position;
        }*/


        //int gridSize = gridSizeNew + (gridSizeNew - 1);

        //var dt = Time.deltaTime;
        //Debug.Log(dt);
        /*
        var mouse = Input.mousePosition;
        var cam = Camera.main;
        var world = cam.ScreenToWorldPoint(new Vector3(mouse.x, mouse.y, cam.nearClipPlane + 30f));
        _particles[Mathf.FloorToInt(gridSize * 0.5f)].position = control.position = world;
        _particles[0].position=(Vector3.zero);*/
    }
    /*void OnGUI()
    {
        //if(Input.GetMouseButtonDown(0))
        if(Event.current.type == EventType.MouseDrag)
        {
            ray = Camera.main.ScreenPointToRay(Input.mousePosition);
            if(Physics.Raycast(ray, out hit))
            {
                var currentScreenPoint = new Vector3(Input.mousePosition.x, Input.mousePosition.y, Input.mousePosition.z);
                var currentPosition = Camera.main.ScreenToWorldPoint(currentScreenPoint);
                var screenPoint = Camera.main.WorldToScreenPoint(transform.position);
                //Debug.Log(currentScreenPoint);
            }
        }

    }*/
    void FindMaxAge(List<Springs> _springs)
    {
        float maxAge = 0f;
        foreach (var p in _springs)
        {
            var dis = p.particleA.Position - p.particleB.Position;
            var dist = dis.magnitude;
            if (dist > maxAge)
            {
                maxAge = dist;
            }
        }
        Debug.Log(maxAge);
    }

    void Update()
    {
        //Chech that a step is done only every deltaTimeStep and press space in the keyboard to pause/start
        timePassed += Time.deltaTime;
        if(timePassed >= deltaTimeStep) timePassed = 0.0f;
        if(!isPaused && timePassed == 0.0f)
        {
            int gridSize = gridSizeNew + gridSizeNew * (gridSizeNew - 1);
            //int gridSize = gridSizeNew + (gridSizeNew - 1);
            simulator.Update(deltaTimeStep, _triangles);
            UpdateMesh();
        }
        if (Input.GetKeyDown(KeyCode.Space))
        {
            isPaused = !isPaused;
        }
        // Sets an anchor with right mouse click
        if (Input.GetMouseButtonDown(1))
        {
            this.ray = Camera.main.ScreenPointToRay(Input.mousePosition);
            if (Physics.Raycast(this.ray, out this.hit))
            {
                this.hit.collider.GetComponent<ParticlesBehaviour>().particles.isActive = false;
            }
            this.screenPoint = Camera.main.WorldToScreenPoint(this.transform.position);
        }

        if(Input.GetKey("left"))
        {
            this.hit.collider.GetComponent<ParticlesBehaviour>().particles.AddPosition(new Vector3(0.05f, 0f, 0f));
        }
        if(Input.GetKey("right"))
        {
            this.hit.collider.GetComponent<ParticlesBehaviour>().particles.AddPosition(new Vector3(-0.05f, 0f, 0f));
        }
        if(Input.GetKey("up"))
        {
            this.hit.collider.GetComponent<ParticlesBehaviour>().particles.AddPosition(new Vector3(0f, 0f, 0.05f));
        }
        if(Input.GetKey("down"))
        {
            this.hit.collider.GetComponent<ParticlesBehaviour>().particles.AddPosition(new Vector3(0f, 0f, -0.05f));
        }
        if(Input.GetKey("w"))
        {
            this.hit.collider.GetComponent<ParticlesBehaviour>().particles.AddPosition(new Vector3(0f, 0.05f, 0f));
        }
        if(Input.GetKey("s"))
        {
            this.hit.collider.GetComponent<ParticlesBehaviour>().particles.AddPosition(new Vector3(0f, -0.05f, 0f));
        }

        // While left mouse click is held down you can drag particles around
        /*if (Input.GetMouseButton(0))
            {
                var currentScreenPoint = new Vector3(Input.mousePosition.x, Input.mousePosition.y, this.screenPoint.z);
                var currentPosition = Camera.main.ScreenToWorldPoint(currentScreenPoint);
                this.hit.collider.GetComponent<ParticlesBehaviour>().particles.Position = currentPosition;
                this.transform.position = currentPosition;
                UpdateMesh();
            }*/

        // Unsets an anchor with middle mouse click
        if (Input.GetMouseButtonDown(2))
        {
            this.ray = Camera.main.ScreenPointToRay(Input.mousePosition);
            if (Physics.Raycast(this.ray, out this.hit))
            {
                this.hit.collider.GetComponent<ParticlesBehaviour>().particles.isActive = true;
            }
        }
        //UpdateMesh();


        /*var mouseRay = Camera.main.ScreenPointToRay(Input.mousePosition);

        RaycastHit hitCollider;
        if (Physics.Raycast(mouseRay, out hitCollider))
        {
            var particlesBehaviour = hitCollider.collider.gameObject.GetComponent<ParticlesBehaviour>();
            if (particlesBehaviour)
            {
                if (Input.GetMouseButtonDown(0))
                {
                    particlesBehaviour.moveWithMouse = true;

                    particlesBehaviour.wasKinematic = particlesBehaviour.particles.isActive;
                    particlesBehaviour.particles.isActive = true;
                }
                if (!Input.GetMouseButton(0) && Input.GetKeyDown(KeyCode.P))
                    particlesBehaviour.particles.isActive = !particlesBehaviour.particles.isActive;
            }
        }*/
    }

    void CreateShape ()
    {
        //int gridSize = gridSizeNew + (gridSizeNew - 1);
        int gridSize = gridSizeNew + gridSizeNew * (gridSizeNew - 1);
        int triangleNum = (gridSize - 1) * (gridSize - 1) * 2 * 3;

        //_triangles = new Triangle[triangleNum];

        int vertexNum = gridSize * gridSize;
        //_vertices = new Vertex[vertexNum];

        triangles = new int[triangleNum];
        vertices = new Vector3 [vertexNum];

        int x = 0;
        for(int j = 0; j < gridSize - 1; j++)
        {
            for(int i = 0; i < gridSize - 1; i++)
            {
                int i0 = j * gridSize + i;
                int i1 = j * gridSize + i + 1;
                int i2 = (j + 1) * gridSize + i;
                int i3 = (j + 1) * gridSize + i +1;


                var p1 = new Triangles(i0, i2, i1);
                p1.PosTriangles(_particles[i0].Position, _particles[i2].Position, _particles[i1].Position);
                _triangles.Add(p1);
                var p2 = new Triangles(i2, i3, i1);
                p2.PosTriangles(_particles[i2].Position, _particles[i3].Position, _particles[i1].Position);
                _triangles.Add(p2);
                triangles[x] = i0;
                x++;
                triangles[x] = i2;
                x++;
                triangles[x] = i1;
                x++;
                triangles[x] = i2;
                x++;
                triangles[x] = i3;
                x++;
                triangles[x] = i1;
                x++;
            }
        } 
        //Debug.Log(_triangles[0].normTri);
        for(int i = 0; i < gridSize; i++)
        {
            for(int j = 0; j < gridSize; j++)
            {
                var pos = j * gridSize + i;
                vertices[pos] = _particles[pos].Position;
            }
        }

    }

    void UpdateMesh()
    {
        //int gridSize = gridSizeNew + (gridSizeNew - 1);
        int gridSize = gridSizeNew + gridSizeNew * (gridSizeNew - 1);
        mesh.Clear();
        int posi = 0;
        for(int j = 0; j < gridSize - 1; j++)
        {
            for(int i = 0; i < gridSize - 1; i++)
            {
                int i0 = j * gridSize + i;
                int i1 = j * gridSize + i + 1;
                int i2 = (j + 1) * gridSize + i;
                int i3 = (j + 1) * gridSize + i +1;

                _triangles[posi].PosTriangles(_particles[i0].Position, _particles[i2].Position, _particles[i1].Position);
                posi++;
                _triangles[posi].PosTriangles(_particles[i2].Position, _particles[i3].Position, _particles[i1].Position);
                posi++;
            }
        } 

        for(int i = 0; i < gridSize; i++)
        {
            for(int j = 0; j < gridSize; j++)
            {
                var pos = j * gridSize + i;
                vertices[pos] = _particles[pos].Position;
            }
        }
        
        mesh.vertices = vertices;
        //mesh.triangles = triangles;

        mesh.SetTriangles(triangles, 0);

        mesh.RecalculateNormals();
        //Debug.Log("Submeshes: " + mesh.subMeshCount);

    }

    void OnDrawGizmos()
    {
        simulator.DrawGizmos();
    }
}