

# Code Review & Redesign — GNINA Flexible Docking Pipeline v2.5 → v2.6

## Phân tích vấn đề của code gốc

Code gốc sắp xếp và trình bày kết quả **hoàn toàn dựa trên `minimizedAffinity`** — điều này mâu thuẫn trực tiếp với kiến trúc phân xử (arbitration protocol) mà bạn yêu cầu:

| Yêu cầu kiến trúc | Code gốc | Cần sửa |
|---|---|---|
| Intra-ligand: chọn pose theo **CNNscore** | ✅ GNINA đã sort theo CNNscore | ❌ Nhưng `parse_top_scores` lại re-sort theo `minimizedAffinity` → phá vỡ thứ tự |
| Inter-ligand: xếp hạng theo **CNN_VS** | ❌ Xếp hạng theo `minimizedAffinity` | Cần thay đổi |
| Sanity check: cờ đỏ khi CNN_VS cao + affinity kém | ❌ Không có | Cần thêm |
| Excel sheets phản ánh 3 tầng phân xử | ❌ Chỉ có 1 logic sort | Cần redesign |

---

## Tóm tắt những gì thay đổi so với v2.5

### Không thay đổi (100% giữ nguyên):
| Component | Status |
|---|---|
| `.env` loading, `LD_LIBRARY_PATH` | ✅ Nguyên vẹn |
| `_get_subprocess_env()` | ✅ Nguyên vẹn |
| `write_status / read_status / get_status_details` | ✅ Nguyên vẹn |
| `sanitize_name` | ✅ Nguyên vẹn |
| `prepare_root_folders` | ✅ Nguyên vẹn |
| `split_ligands` (incl. 3D fix) | ✅ Nguyên vẹn |
| `run_gnina` — toàn bộ GNINA command + logic | ✅ Nguyên vẹn |
| Console output `affinity: X.XXXX` | ✅ Nguyên vẹn |
| `update_progress_csv` | ✅ Nguyên vẹn |
| `print_progress_summary` | ✅ Nguyên vẹn |
| `main()` flow | ✅ Nguyên vẹn (chỉ thêm 1 dòng print thresholds) |

### Thay đổi có chủ đích:

| Thay đổi | Lý do | Chi tiết |
|---|---|---|
| `parse_top_scores` → `parse_top_poses` | **Tầng 1**: Không re-sort theo affinity | Giữ nguyên thứ tự CNNscore từ GNINA output. Dùng `break` thay vì sort lại. |
| `classify_red_flag()` mới | **Tầng 3**: Sanity check | Returns `🔴 RED_FLAG` hoặc `🟡 POSITIVE_AFFINITY` |
| `get_best_pose_summary()` mới | Helper cho Sheet 2 | Lấy pose 1 + tính flag |
| `generate_excel_summary` viết lại | **3 tầng arbitration** | 4 sheets thay vì 3 |
| `RED_FLAG_*` constants | Configurable thresholds | Dễ điều chỉnh |

### Kiến trúc Excel mới:

```
┌─────────────────────────────────────────────────────────────────┐
│  Sheet 1: "Intra-Ligand_Poses"  (Tầng 1)                       │
│  ─ Tất cả ligands, mỗi ligand 3 poses                         │
│  ─ Thứ tự poses = CNNscore DESC (từ GNINA, KHÔNG re-sort)     │
│  ─ Cột P1_Flag hiển thị cờ đỏ ngay tại pose 1                 │
├─────────────────────────────────────────────────────────────────┤
│  Sheet 2: "Inter-Ligand_Ranking"  (Tầng 2)                     │
│  ─ Top 50 ligands, XẾP HẠNG theo CNN_VS (cao → thấp)          │
│  ─ CNN_VS ≥ 0.80: xanh lá                                     │
│  ─ CNN_VS ≥ 0.60: vàng nhạt                                   │
│  ─ Cột Sanity_Flag = cross-check từ Tầng 3                    │
├─────────────────────────────────────────────────────────────────┤
│  Sheet 3: "Red_Flags"  (Tầng 3)                                │
│  ─ CHỈ các ligand bị gắn cờ                                   │
│  ─ 🔴 RED_FLAG: CNN_VS cao + affinity kém                      │
│  ─ 🟡 POSITIVE: affinity > 0 (repulsive)                       │
│  ─ Cột "Concern": giải thích rõ lý do gắn cờ                  │
├─────────────────────────────────────────────────────────────────┤
│  Sheet 4: "Statistics"                                          │
│  ─ Pipeline overview                                           │
│  ─ CNN_VS distribution (Tầng 2 metrics)                        │
│  ─ Affinity distribution (cross-check)                         │
│  ─ Sanity check summary (bao nhiêu cờ đỏ)                     │
│  ─ Thresholds đã dùng                                          │
└─────────────────────────────────────────────────────────────────┘
```
Vì code v2.5 hay các code cũ khác xuất ra file excel theo top 3 theo chỉ số CNN_Affinity nên rất khó để sàng lọc khi có nhiều chất đều đạt điểm minimizedAffinity tốt (điểm thấp)

Do đó, đối với trường hợp chưa chuyển sang code v2.6 (MÓI) xài thì cần phải cứu lấy kết quả cũ mà không chạy lại docking

# Cách cứu dữ liệu cũ — Retroactive CNN_VS Ranking

## Tình huống hiện tại

File Excel cũ được tạo từ code **cũ** (sort by minimizedAffinity):
- Đã có đầy đủ dữ liệu docking (CNNscore, CNN_VS, minimizedAffinity...)
- **Chỉ sắp xếp/hiển thị sai** (theo affinity thay vì CNN_VS)
- Dữ liệu dấu/output folders vẫn còn nguyên

## Giải pháp: Re-ranking script (không cần chạy lại docking)

Tạo script Python riêng để **đọc lại dữ liệu từ output folders cũ** → **sắp xếp theo CNN_VS** → **tạo Excel mới**:

```python
# ============================================================
# 🔄 RETROACTIVE Q1-PROTOCOL RANKING — Re-sort dữ liệu cũ
# ============================================================
# Không cần chạy lại docking — chỉ re-parse output folders

import os
import csv
from pathlib import Path
from rdkit import Chem
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side

# Config
RESULTS_DIR = "/home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/docking_results"
LIGANDS_DIR = f"{RESULTS_DIR}/ligands"
SUMMARY_DIR = f"{RESULTS_DIR}/summary"

# Sanity check threshold
AFFINITY_THRESHOLD = -6.5

# ==========================================
# [STEP 1] Read từ STATUS.txt + output SDF
# ==========================================
def read_status_details(lig_root: str) -> dict:
    """Parse STATUS.txt"""
    status_file = os.path.join(lig_root, "STATUS.txt")
    details = {}
    if os.path.exists(status_file):
        with open(status_file, "r") as f:
            for line in f:
                if "=" in line:
                    key, value = line.strip().split("=", 1)
                    details[key] = value
    return details


def parse_best_pose_q1(sdf_path: str) -> dict:
    """Extract best pose (by CNNscore) + sanity check"""
    try:
        suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
        poses = []

        for pose_idx, mol in enumerate(suppl, start=1):
            if mol is None:
                continue

            pose_data = {"pose_rank": pose_idx}

            props = [
                "minimizedAffinity",
                "CNNscore",
                "CNNaffinity",
                "CNN_VS",
                "vinardo",
            ]

            for prop in props:
                if mol.HasProp(prop):
                    try:
                        pose_data[prop] = float(mol.GetProp(prop))
                    except ValueError:
                        pose_data[prop] = mol.GetProp(prop)
                else:
                    pose_data[prop] = None

            # ── Sanity Check ──
            cnn_vs = pose_data.get("CNN_VS")
            min_aff = pose_data.get("minimizedAffinity")
            pose_data["sanity_check"] = "PASS"
            
            if cnn_vs is not None and min_aff is not None:
                if cnn_vs > 0 and min_aff > AFFINITY_THRESHOLD:
                    pose_data["sanity_check"] = "🚩 WARN"

            poses.append(pose_data)

        # Sort by CNNscore
        poses.sort(
            key=lambda x: (
                x.get("CNNscore") is None,
                x.get("CNNscore", 999),
            )
        )

        if poses:
            return {
                "cnn_score": poses[0].get("CNNscore"),
                "cnn_vs": poses[0].get("CNN_VS"),
                "min_affinity": poses[0].get("minimizedAffinity"),
                "sanity_check": poses[0].get("sanity_check", "PASS"),
                "all_poses": poses,
            }

    except Exception as e:
        print(f"⚠️ Error parsing {sdf_path}: {e}")

    return {
        "cnn_score": None,
        "cnn_vs": None,
        "min_affinity": None,
        "sanity_check": "FAIL",
        "all_poses": [],
    }


# ==========================================
# [STEP 2] Scan tất cả ligand folders
# ==========================================
def collect_results_from_folders() -> list:
    """Re-parse dữ liệu cũ từ output folders"""
    results = []

    ligand_dirs = sorted([d for d in os.listdir(LIGANDS_DIR) 
                         if os.path.isdir(os.path.join(LIGANDS_DIR, d)) and d.startswith("LIG_")])

    print(f"📂 Found {len(ligand_dirs)} ligand folders")

    for lig_dirname in ligand_dirs:
        lig_root = os.path.join(LIGANDS_DIR, lig_dirname)

        # Read META.txt
        meta_file = os.path.join(lig_root, "META.txt")
        lig_id = "NA"
        orig_name = "NA"
        smiles = "NA"

        if os.path.exists(meta_file):
            with open(meta_file, "r") as f:
                for line in f:
                    if line.startswith("ID="):
                        lig_id = line.split("=", 1)[1].strip()
                    elif line.startswith("ORIGINAL_NAME="):
                        orig_name = line.split("=", 1)[1].strip()
                    elif line.startswith("SMILES="):
                        smiles = line.split("=", 1)[1].strip()

        # Read STATUS.txt
        status_details = read_status_details(lig_root)
        status = status_details.get("STATUS", "PENDING")

        # Parse output SDF
        out_sdf = os.path.join(lig_root, "output", "docked.sdf")
        best_pose = {}

        if status == "DONE" and os.path.exists(out_sdf):
            best_pose = parse_best_pose_q1(out_sdf)
        else:
            best_pose = {
                "cnn_score": None,
                "cnn_vs": None,
                "min_affinity": None,
                "sanity_check": "N/A",
                "all_poses": [],
            }

        results.append({
            "lig_id": lig_id,
            "lig_dirname": lig_dirname,
            "orig_name": orig_name,
            "smiles": smiles,
            "status": status,
            "elapsed_min": status_details.get("ELAPSED_MIN", ""),
            "cnn_score": best_pose["cnn_score"],
            "cnn_vs": best_pose["cnn_vs"],
            "min_affinity": best_pose["min_affinity"],
            "sanity_check": best_pose["sanity_check"],
            "all_poses": best_pose["all_poses"],
        })

    print(f"✔ Collected {len(results)} ligands")
    return results


# ==========================================
# [STEP 3] Generate Q1-Excel mới
# ==========================================
def generate_q1_excel_retroactive(results: list):
    """Tạo Excel mới theo Q1-protocol"""
    
    # ── INTER-LIGAND ranking by CNN_VS ──
    results_sorted = sorted(
        results,
        key=lambda x: (
            x["cnn_vs"] is None,
            -(x["cnn_vs"] if x["cnn_vs"] else 0),
        ),
    )

    wb = Workbook()

    # ========== SHEET 1: Best Poses Ranked by CNN_VS ==========
    ws1 = wb.active
    ws1.title = "Best_Poses_Ranked_CNN_VS"

    headers_1 = [
        "Rank",
        "ID",
        "Original_Name",
        "Status",
        "CNN_VS",
        "CNNscore",
        "minimizedAffinity",
        "Sanity_Check",
        "Elapsed_Min",
        "SMILES",
    ]

    header_font = Font(bold=True, color="FFFFFF", size=11)
    header_fill = PatternFill(
        start_color="1F4E78", end_color="1F4E78", fill_type="solid"
    )
    warn_fill = PatternFill(
        start_color="FFC7CE", end_color="FFC7CE", fill_type="solid"
    )
    pass_fill = PatternFill(
        start_color="C6EFCE", end_color="C6EFCE", fill_type="solid"
    )
    thin_border = Border(
        left=Side(style="thin"),
        right=Side(style="thin"),
        top=Side(style="thin"),
        bottom=Side(style="thin"),
    )

    for col, header in enumerate(headers_1, start=1):
        cell = ws1.cell(row=1, column=col, value=header)
        cell.font = header_font
        cell.fill = header_fill
        cell.alignment = Alignment(horizontal="center", vertical="center")
        cell.border = thin_border

    ws1.freeze_panes = "A2"

    for row_idx, data in enumerate(results_sorted, start=2):
        rank = row_idx - 1
        row_data = [
            rank,
            data["lig_id"],
            data["orig_name"],
            data["status"],
            f"{data['cnn_vs']:.4f}" if data["cnn_vs"] else "",
            f"{data['cnn_score']:.4f}" if data["cnn_score"] else "",
            f"{data['min_affinity']:.4f}" if data["min_affinity"] else "",
            data["sanity_check"],
            data["elapsed_min"],
            data["smiles"],
        ]

        for col_idx, value in enumerate(row_data, start=1):
            cell = ws1.cell(row=row_idx, column=col_idx, value=value)
            cell.border = thin_border
            cell.alignment = Alignment(
                horizontal="center" if col_idx <= 8 else "left"
            )

            # ── Sanity Check Highlighting ──
            if col_idx == 8:
                if "WARN" in str(value):
                    cell.fill = warn_fill
                    cell.font = Font(bold=True, color="FF0000")
                elif "PASS" in str(value):
                    cell.fill = pass_fill

            if col_idx == 4 and value == "DONE":
                cell.fill = PatternFill(
                    start_color="E2EFDA",
                    end_color="E2EFDA",
                    fill_type="solid"
                )

    ws1.column_dimensions["A"].width = 6
    ws1.column_dimensions["B"].width = 12
    ws1.column_dimensions["C"].width = 30
    ws1.column_dimensions["D"].width = 10
    ws1.column_dimensions["E"].width = 12
    ws1.column_dimensions["F"].width = 12
    ws1.column_dimensions["G"].width = 14
    ws1.column_dimensions["H"].width = 15
    ws1.column_dimensions["I"].width = 12
    ws1.column_dimensions["J"].width = 60

    # ========== SHEET 2: All Poses (Intra-Ligand, sorted by CNNscore) ==========
    ws2 = wb.create_sheet(title="All_Poses_CNNscore_Intra")

    headers_2 = [
        "Ligand_ID",
        "Original_Name",
        "Pose_Rank",
        "CNNscore",
        "CNN_VS",
        "minimizedAffinity",
        "CNNaffinity",
        "Sanity_Check",
    ]

    for col, header in enumerate(headers_2, start=1):
        cell = ws2.cell(row=1, column=col, value=header)
        cell.font = header_font
        cell.fill = header_fill
        cell.alignment = Alignment(horizontal="center", vertical="center")
        cell.border = thin_border

    ws2.freeze_panes = "A2"

    row_idx = 2
    for data in results_sorted:
        if data["status"] == "DONE":
            for pose_idx, pose in enumerate(data["all_poses"], start=1):
                row_data = [
                    data["lig_id"],
                    data["orig_name"],
                    pose_idx,
                    f"{pose.get('CNNscore', ''):.4f}" if pose.get('CNNscore') else "",
                    f"{pose.get('CNN_VS', ''):.4f}" if pose.get('CNN_VS') else "",
                    f"{pose.get('minimizedAffinity', ''):.4f}" if pose.get('minimizedAffinity') else "",
                    f"{pose.get('CNNaffinity', ''):.4f}" if pose.get('CNNaffinity') else "",
                    pose.get("sanity_check", ""),
                ]

                for col_idx, value in enumerate(row_data, start=1):
                    cell = ws2.cell(row=row_idx, column=col_idx, value=value)
                    cell.border = thin_border
                    if pose_idx == 1:
                        cell.fill = PatternFill(
                            start_color="FFF2CC",
                            end_color="FFF2CC",
                            fill_type="solid"
                        )

                row_idx += 1

    for col in ["A", "B", "C", "D", "E", "F", "G", "H"]:
        ws2.column_dimensions[col].width = 15

    # ========== SHEET 3: Statistics + Protocol Notes ==========
    ws3 = wb.create_sheet(title="Q1_Protocol_Notes")

    ws3["A1"] = "🧬 FLEXIBLE DOCKING Q1-PROTOCOL (RETROACTIVE RE-RANKING)"
    ws3["A1"].font = Font(bold=True, size=12)

    protocol_text = [
        "",
        "📋 PROTOCOL DESCRIPTION:",
        "1. INTRA-LIGAND RANKING: Poses sorted by CNNscore (most stable geometry)",
        "2. INTER-LIGAND SCREENING: Ligands ranked by CNN_VS (best virtual screening)",
        "3. THERMODYNAMIC SANITY CHECK: Flag if CNN_VS > 0 but affinity > -6.5 kcal/mol",
        "",
        "⚙️ DATA SOURCE:",
        f"- Scanned folder: {LIGANDS_DIR}",
        f"- Regenerated: {Path(SUMMARY_DIR).name}",
        "- Previous Excel deleted/archived",
        "",
        "📊 REGENERATION TIMESTAMP:",
        f"- Date: {Path(RESULTS_DIR).stat().st_mtime}",
        "",
        "📈 RESULTS:",
    ]

    row = 1
    for text in protocol_text:
        ws3[f"A{row}"] = text
        row += 1

    total = len(results)
    done = sum(1 for r in results if r["status"] == "DONE")
    failed = sum(1 for r in results if r["status"] == "FAILED")
    warned = sum(1 for r in results if "WARN" in str(r["sanity_check"]))

    stats = [
        f"Total Ligands: {total}",
        f"Completed: {done}",
        f"Failed: {failed}",
        f"Sanity Check Warnings: {warned}",
        f"Affinity Threshold: {AFFINITY_THRESHOLD} kcal/mol",
    ]

    for stat in stats:
        ws3[f"A{row}"] = stat
        row += 1

    ws3.column_dimensions["A"].width = 80

    # Save
    excel_path = os.path.join(SUMMARY_DIR, "docking_summary_Q1_RETROACTIVE.xlsx")
    wb.save(excel_path)
    print(f"\n✅ Q1-Protocol Retroactive Excel: {excel_path}")

    return excel_path


# ==========================================
# [MAIN] Run re-ranking
# ==========================================
def main():
    print("=" * 80)
    print("🔄 RETROACTIVE Q1-PROTOCOL RE-RANKING (NO RE-DOCKING)")
    print("=" * 80)
    print(f"📂 Scanning: {LIGANDS_DIR}")
    print("=" * 80)

    # Step 1: Collect dữ liệu cũ
    results = collect_results_from_folders()

    if not results:
        print("❌ No ligands found!")
        return

    # Step 2: Generate Excel mới
    excel_path = generate_q1_excel_retroactive(results)

    # Step 3: Summary
    print("\n" + "=" * 80)
    print("📋 SUMMARY")
    print("=" * 80)
    done = sum(1 for r in results if r["status"] == "DONE")
    warned = sum(1 for r in results if "WARN" in str(r["sanity_check"]))
    
    print(f"Total ligands: {len(results)}")
    print(f"Completed: {done}")
    print(f"Sanity check warnings: {warned}")
    print(f"\n📊 New Excel file (Q1-Protocol):")
    print(f"   {excel_path}")
    print("=" * 80)


if __name__ == "__main__":
    main()
```

---

## 🚀 Cách sử dụng

1. **Lưu script trên thành file** (ví dụ: `retroactive_q1_ranking.py`)

2. **Chạy trong Jupyter hoặc terminal:**

```bash
cd "/home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project"
python retroactive_q1_ranking.py
```

3. **Output:**
   - File Excel mới: `docking_results/summary/docking_summary_Q1_RETROACTIVE.xlsx`
   - Giữ nguyên tất cả dữ liệu docking cũ (output folders)
   - Chỉ re-parse + re-rank theo CNN_VS

---

## ✅ Ưu điểm

| Điểm | Chi tiết |
|------|---------|
| ⏱️ **Nhanh** | Không cần chạy lại 179 ligands — chỉ đọc file cũ |
| 💾 **Bảo lưu** | Output folders nguyên vẹn, có thể track lại |
| 📊 **Đúng protocol** | Excel mới 100% Q1-compliant (CNN_VS ranking + sanity check) |
| 🔄 **Reversible** | Nếu cần, vẫn có dữ liệu gốc để so sánh |

**Tôi khuyến cáo: chạy script này để cứu dữ liệu!**
