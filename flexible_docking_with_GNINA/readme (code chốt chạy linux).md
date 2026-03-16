

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
