import { useEffect, useRef } from "react";
import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { buildStructureGeometry, type Vec3 } from "../lib/crystal";
import type { Material } from "../data/materials";

interface Props {
  material: Material;
}

function centerOf(points: Vec3[]): Vec3 {
  const c: Vec3 = [0, 0, 0];
  for (const p of points) { c[0] += p[0]; c[1] += p[1]; c[2] += p[2]; }
  return [c[0] / points.length, c[1] / points.length, c[2] / points.length];
}

export default function StructureViewer({ material }: Props) {
  const mountRef = useRef<HTMLDivElement>(null);
  const stateRef = useRef<{
    renderer: THREE.WebGLRenderer;
    scene: THREE.Scene;
    camera: THREE.PerspectiveCamera;
    controls: OrbitControls;
    group: THREE.Group;
    raf: number;
  } | null>(null);

  // one-time scene setup
  useEffect(() => {
    const mount = mountRef.current;
    if (!mount) return;

    const scene = new THREE.Scene();
    scene.background = null;

    const camera = new THREE.PerspectiveCamera(40, 1, 0.1, 100);
    camera.position.set(2.1, 1.7, 2.4);

    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    mount.appendChild(renderer.domElement);

    const controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.08;
    controls.autoRotate = true;
    controls.autoRotateSpeed = 1.6;
    controls.minDistance = 1.2;
    controls.maxDistance = 8;

    scene.add(new THREE.AmbientLight(0xffffff, 0.75));
    const key = new THREE.DirectionalLight(0xffffff, 1.1);
    key.position.set(3, 4, 2);
    scene.add(key);
    const fill = new THREE.DirectionalLight(0xffffff, 0.4);
    fill.position.set(-3, -2, -3);
    scene.add(fill);

    const group = new THREE.Group();
    scene.add(group);

    const resize = () => {
      const w = mount.clientWidth;
      const h = mount.clientHeight;
      if (w === 0 || h === 0) return;
      camera.aspect = w / h;
      camera.updateProjectionMatrix();
      renderer.setSize(w, h);
    };
    resize();
    const ro = new ResizeObserver(resize);
    ro.observe(mount);

    let raf = 0;
    const animate = () => {
      controls.update();
      renderer.render(scene, camera);
      raf = requestAnimationFrame(animate);
    };
    animate();

    stateRef.current = { renderer, scene, camera, controls, group, raf };

    return () => {
      cancelAnimationFrame(raf);
      ro.disconnect();
      controls.dispose();
      renderer.dispose();
      mount.removeChild(renderer.domElement);
      stateRef.current = null;
    };
  }, []);

  // rebuild geometry whenever the material changes
  useEffect(() => {
    const st = stateRef.current;
    if (!st) return;
    const { group } = st;
    while (group.children.length) {
      const child = group.children[0];
      group.remove(child);
      if ((child as THREE.Mesh).geometry) (child as THREE.Mesh).geometry.dispose();
    }

    const geo = buildStructureGeometry(material);
    const allPoints = geo.cellEdges.flat();
    const center = centerOf(allPoints);
    const toV3 = (p: Vec3) => new THREE.Vector3((p[0] - center[0]) * geo.scale, (p[1] - center[1]) * geo.scale, (p[2] - center[2]) * geo.scale);

    const sphereGeo = new THREE.SphereGeometry(1, 20, 16);
    for (const atom of geo.atoms) {
      const radius = atom.role === "cation" ? 0.085 : 0.065;
      const mat = new THREE.MeshStandardMaterial({ color: atom.color, roughness: 0.45, metalness: 0.15 });
      const mesh = new THREE.Mesh(sphereGeo, mat);
      mesh.scale.setScalar(radius);
      const p = toV3(atom.position);
      mesh.position.copy(p);
      group.add(mesh);
    }

    const bondMat = new THREE.LineBasicMaterial({ color: 0x9aa39a, transparent: true, opacity: 0.6 });
    const bondPositions: number[] = [];
    for (const [a, b] of geo.bonds) {
      const pa = toV3(a);
      const pb = toV3(b);
      bondPositions.push(pa.x, pa.y, pa.z, pb.x, pb.y, pb.z);
    }
    const bondGeo = new THREE.BufferGeometry();
    bondGeo.setAttribute("position", new THREE.Float32BufferAttribute(bondPositions, 3));
    group.add(new THREE.LineSegments(bondGeo, bondMat));

    const edgeMat = new THREE.LineBasicMaterial({ color: 0x9aa39a, transparent: true, opacity: 0.28 });
    const edgePositions: number[] = [];
    for (const [a, b] of geo.cellEdges) {
      const pa = toV3(a);
      const pb = toV3(b);
      edgePositions.push(pa.x, pa.y, pa.z, pb.x, pb.y, pb.z);
    }
    const edgeGeo = new THREE.BufferGeometry();
    edgeGeo.setAttribute("position", new THREE.Float32BufferAttribute(edgePositions, 3));
    group.add(new THREE.LineSegments(edgeGeo, edgeMat));
  }, [material]);

  return <div ref={mountRef} className="structure-canvas" />;
}
