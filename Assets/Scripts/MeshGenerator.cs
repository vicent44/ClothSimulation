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

    [SerializeField] float elasticConstant = 1f;
    [SerializeField] float dampingConstant = 1f;

    Simulate simulator;
    TriangleIntersection intersec;
    Hashing hash;

    //Triangles triangle;
    List<Particles> _particles;
    List<Springs> _springs;
    List<Triangles> _triangles;
    //For hash table

    public struct SPHash
    {
        public List<int> indices;
    } 

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

        //float mass = massTotalCloth / (gridSize * gridSize);
        //float area = areaToalCloth / (gridSize * gridSize);
        float mass = 0.1f;

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
                var p = new Particles(new Vector3(jP,0.0f,iP), mass, pi, pj);
                _particles.Add(p);
                //var obj = GameObject.Instantiate(_esferePrefab, new Vector3(jP,0.0f,iP), Quaternion.identity);
                var obj = GameObject.CreatePrimitive(PrimitiveType.Sphere);
                obj.transform.localScale = Vector3.one * 0.1f;
                obj.transform.position = p.Position;
                //obj.AddComponent<Particles>();
                //obj.gameObject.tag = "("+ jP.ToString() +" ,"+ y.ToString() +" ,"+ iP.ToString() +")";
                _esfere.Add(obj);
                pi++;
            }
            pj++;
        }
        float distBetwenParticles = (_particles[0].Position - _particles[1].Position).magnitude;
        float area = (distBetwenParticles * distBetwenParticles * (gridSize - 1) * (gridSize - 1)) / (gridSize * gridSize);
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
                var re = new Springs(c, right, elasticConstant, dampingConstant, 1);
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

                var re = new Springs(c, right, elasticConstant, dampingConstant, 1);
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

                var re = new Springs(c, bot, elasticConstant, dampingConstant, 2);
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
                var re = new Springs(c, top, elasticConstant, dampingConstant, 2);
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

                var re = new Springs(c, left, elasticConstant, dampingConstant, 3);
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

                var re = new Springs(c, left, elasticConstant, dampingConstant, 3);
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
        _particles[15].isActive = false;

        simulator = new Simulate(_particles, _springs, _triangles, windforce, plane);
        intersec = new TriangleIntersection();

        mesh = new Mesh();
        GetComponent<MeshFilter>().mesh = mesh;
        CreateShape();
        UpdateMesh();


    }

    void FixedUpdate()
    {
        int gridSize = gridSizeNew + gridSizeNew * (gridSizeNew - 1);
        simulator.Update(Time.fixedDeltaTime, _triangles);

        //Collision constrain

        hash = new Hashing(gridSize, 1.0f/(float)gridSize, 2000);
        SPHash[] spHash = new SPHash[2000];

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
                        if(idx != _particles[tri.indexTriA] && idx != _particles[tri.indexTriB] && idx != _particles[tri.indexTriC])
                        {
                            Vector3 p = _particles[idx].Position;
                            float w = _particles[idx].Mass;

                            Vector3 corr, corr0, corr1, corr2;
                            if(TrianglePointDistanceConstraint(
                                p, w,
                                p0, w0,
                                p1, w1,
                                p2, w2,
                                thickness, 1f, 0.0f,
                                out corr, out corr0, out corr1, out corr2))
                            {
                                _particles[idx].Position += corr;
                                _particles[tri.indexTriA].Position += corr0;
                                _particles[tri.indexTriB].Position += corr1;
                                _particles[tri.indexTriC].Position += corr2;
                            }
                        }
                    }
                }
            }
        }

        UpdateMesh();

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


        for(int i = 0; i < _particles.Count; i++)
        {
            _esfere[i].transform.position = _particles[i].Position;
        }


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