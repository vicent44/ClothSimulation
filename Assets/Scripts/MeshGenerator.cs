﻿using System.Collections;
using System.Collections.Generic;
using UnityEngine;


[RequireComponent(typeof(MeshFilter))]
public class MeshGenerator : MonoBehaviour
{
    Mesh mesh;

    Vector3[] vertices;
    Vector3[] verticesNext;
    int[] triangles;
    public int gridSizeNew;
    //public bool gravity = true;
    [SerializeField] int iterations = 1;

    //[SerializeField] float areaTotalCloth = 12f;

    [SerializeField] float massTotalCloth = 12f;
    [SerializeField] float clothDensity = 1f;


    //[SerializeField] bool useGravity = true;
    //[SerializeField] float gravity = 1f;

    //[SerializeField] bool useWind = false;
    //[SerializeField] Vector3 Vector3Field (Rect position, string windDirection, Vector3(0, 0, 1));
    /*[SerializeField] int windxDirection = 1;
    [SerializeField] int windyDirection = 0;
    [SerializeField] int windzDirection = 0;*/
    [SerializeField] Vector3 windDirection=new Vector3(0,0,1);

    [SerializeField] float windModule = 1f;

    //[SerializeField] float elasticConstant = 1f;
    [SerializeField] float elasticConstantStructural = 20f;
    [SerializeField] float elasticConstantShear = 20f;
    [SerializeField] float elasticConstantBend = 20f;
    [SerializeField] float dampingConstant = 0.25f;

    Simulate simulator;
    TriangleIntersection intersec;
    Hashing hash;

    float timePassed = 0.0f;
    public float deltaTimeStep = 0.02f;

    //Triangles triangle;
    List<Particles> _particles;
    List<Springs> _springs;
    List<Triangles> _triangles;
    //For hash table
/*
    public struct SPHash
    {
        public List<int> indices;
    } */

    Vector3 windforce;

    //[SerializeField] Transform control;
    //[SerializeField] Transform _esferePrefab;
    //[SerializeField] GameObject esferePrefab;

    [SerializeField] GameObject esfereControl;
    List<GameObject> _esfere;
    //List<Transform> _esfere;

    [SerializeField] Transform plane;

    //Structural = 1
    //Shear = 2
    //Bend = 3

    //Mouse Drag Variables
    Vector3 screenPoint;
    Vector3 offset;
    Ray ray;
    RaycastHit hit;

    public bool isPaused = true;

    void Start()
    {
        /*mesh = new Mesh();
        GetComponent<MeshFilter>().mesh = mesh;
        CreateShape();
        UpdateMesh();*/
        //int gridSize = gridSizeNew + (gridSizeNew - 1);
        
        int gridSize = gridSizeNew + gridSizeNew * (gridSizeNew - 1);
        //int gridSize = 4;
        _particles = new List<Particles>();
        _springs = new List<Springs>();
        _triangles = new List<Triangles>();
        //_esfere = new List<Transform>();
        _esfere = new List<GameObject>();
        //windforce = new Vector3();


        var particleParent = new GameObject("particleParent");

        //var particle = new GameObject();

        float area = 0.2f * (gridSize - 1) * 0.2f * (gridSize - 1);

        //float mass = massTotalCloth / (gridSize * gridSize);
        //float area = areaToalCloth / (gridSize * gridSize);
        //float mass = 0.1f;
        float massTotal = clothDensity * area;
        float mass = massTotal / (gridSize * gridSize);

        //int pos = 0;
        float jP = 0;
        float iP = 0;
        int pi = 0;
        int pj = 0;
        //float y = 0;
        int indexat = 0;
        for(int i = 0; i < gridSize; i++)
        {
            for(int j = 0; j < gridSize; j++)
            {
                //pos = j * gridSize + i;
                jP = j * 0.2f;
                iP = i * 0.2f;
                
                //var p = new Particles(new Vector3(jP,0.0f,iP), mass, pi, pj);
                var newP = ParticlesBehaviour.Create(new Vector3(jP,0.0f,iP), mass, pi, pj, esfereControl);
                newP.transform.SetParent(particleParent.transform, false);
                newP.transform.localPosition = new Vector3(jP,0.0f,iP);
                newP.particles.Position = newP.transform.position;

                //newP.transform.SetParent(particleParent, false);
                //newP.transform.localPosition = new Vector3(jP,0.0f,iP);
                //_esfere[pi].AddComponent<p>();
                //p.particleObject = Instantiate(_esferePrefab, p.Position, Quaternion.identity);//AddComponent<SphereCollider>();//<SphereCollider>();
                //_esfere.Add(newP);
                _particles.Add(newP.particles);
                //var obj = GameObject.Instantiate(_esferePrefab, new Vector3(jP,0.0f,iP), Quaternion.identity);
                //var obj = GameObject.CreatePrimitive(PrimitiveType.Sphere);
                //obj.transform.localScale = Vector3.one * 0.1f;
                //obj.transform.position = p.Position;
                //obj.AddComponent<Particles>();
                //obj.gameObject.tag = "("+ jP.ToString() +" ,"+ y.ToString() +" ,"+ iP.ToString() +")";
                //_esfere.Add(obj);
                pi++;
            }
            pj++;
        }
        float distBetwenParticles = (_particles[0].Position - _particles[1].Position).magnitude;
        //float area = (distBetwenParticles * distBetwenParticles * (gridSize - 1) * (gridSize - 1)) / (gridSize * gridSize);
        //var w = new WindForce(new Vector3(windxDirection, windyDirection, windzDirection) , windModule, area, clothDensity);
        //windforce = w.WindTotalForce;
        //windforce = new Vector3((float)windxDirection, (float)windyDirection, (float)windzDirection) * windModule * area * clothDensity;
        //windforce.x = (float)windxDirection * windModule * area * clothDensity;
        windforce = windDirection * windModule * area * clothDensity;
        //Debug.Log(windforce.z);
        //Springs

        //Structural Horizontal - Right
        for(int j = 0; j < gridSize; j++)
        {
            for(int i = 0; i < gridSize - 1; i++)
            {
                var index = j * gridSize + i;
                var c = _particles[index];
                var right = _particles[index + 1];
                //Debug.Log(index);
                var re = new Springs(c, right, elasticConstantStructural, dampingConstant, 1);
                //c.Connect(re);
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
                //c.Connect(re);
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
                //c.Connect(re);
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
                //Debug.Log(index_1);
                var re = new Springs(c, top, elasticConstantShear, dampingConstant, 2);
                //c.Connect(re);
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
                //c.Connect(re);
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
                //c.Connect(re);
                _springs.Add(re);
            }
        }
        
        //Fix first row
        /*for(int i = 0; i < 16; i++)
        {
            _particles[i].isActive = false;
        }*/
        _particles[0].isActive = false;
        _particles[8].isActive = false;
        //_particles[3].isActive = false;

        simulator = new Simulate(_particles, _springs, _triangles, windforce, plane, gridSize);
        //intersec = new TriangleIntersection();

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
        timePassed += Time.deltaTime;
        if(timePassed >= deltaTimeStep) timePassed = 0.0f;
        if(!isPaused && timePassed == 0.0f)
        {
            int gridSize = gridSizeNew + gridSizeNew * (gridSizeNew - 1);
            simulator.Update(deltaTimeStep, _triangles);
            UpdateMesh();
            //Debug.Log(deltaTimeStep);
        }
        if (Input.GetKeyDown(KeyCode.Space))
            isPaused = !isPaused;

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
            //Debug.Log("hola");
            //_particles[0].isActive = true;
            this.hit.collider.GetComponent<ParticlesBehaviour>().particles.AddPosition(new Vector3(0.05f, 0f, 0f));
            //_particles[0].ResetResultantForce();
            //_particles[0].Velocity = Vector3.zero;
            //_particles[0].isActive = false;
            
            //_particles[0].ResetResultantForce();
        }
        if(Input.GetKey("right"))
        {
            //_particles[0].isActive = true;
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
        //Debug.Log(_triangles[30].normTri);     

        for(int i = 0; i < gridSize; i++)
        {
            for(int j = 0; j < gridSize; j++)
            {
                var pos = j * gridSize + i;
                vertices[pos] = _particles[pos].Position;
            }
        }
        
        mesh.vertices = vertices;
        mesh.triangles = triangles;

        mesh.RecalculateNormals();

    }

    void OnDrawGizmos()
    {
        simulator.DrawGizmos();
    }
}