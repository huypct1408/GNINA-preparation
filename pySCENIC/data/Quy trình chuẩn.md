**LỆNH ĐÓNG BĂNG (FREEZE COMMAND) ĐÃ ĐƯỢC KÍCH HOẠT.** Toàn bộ khung lý thuyết đã đạt đến trạng thái bão hòa hoàn hảo. Không thêm bất kỳ thuật toán hay biến số nào nữa. 

Dưới đây là **Quy trình Thực thi Tiêu chuẩn (Standard Operating Procedure - SOP)** duy nhất và cuối cùng cho Kiến trúc Dược lý Hệ thống 4 Lớp của bạn. Quy trình này được thiết kế theo chuẩn Q1, kết hợp chặt chẽ hệ quy chiếu **SMART**, tiêu chuẩn quản trị dữ liệu **FAIR**, và đặc biệt là bộ **"Kỷ luật Thép" (Red Flags)** đúc kết từ 10 lỗ hổng phản biện nhằm ngăn chặn mọi sự sụp đổ về mặt nhân quả.

---

### 🚫 CÁC TÌNH HUỐNG TỐI KỴ BẮT BUỘC PHẢI TRÁNH (FATAL ANTI-PATTERNS)
Nếu vi phạm một trong các điều sau, toàn bộ pipeline sẽ xuất ra "Rác" (Garbage In, Garbage Out) và bị từ chối xuất bản:
1. **Tuyệt đối không chạy RWR mà không qua Cổng Nhiệt động học:** Không đưa các hợp chất có $\Delta G \ge -6.5$ kcal/mol vào Lớp 2. Việc chuẩn hóa $P_0$ trên các chất yếu sẽ gây ra "Nghịch lý Chuẩn hóa", làm thuật toán RWR nhầm lẫn chất vô tác dụng với siêu thuốc.
2. **Tuyệt đối không dùng Ma trận Vô hướng (Undirected Graph) cho toàn mạng lưới:** Mạng lưới phiên mã (pySCENIC) BẮT BUỘC phải là ma trận có hướng (`nx.DiGraph`). Cho phép tín hiệu chảy ngược từ gen lên Yếu tố phiên mã (TF) sẽ vi phạm Nguyên lý Trung tâm của Sinh học.
3. **Tuyệt đối không gọi HEK-293 là "Tế bào khỏe mạnh bình thường":** HEK-293 là tế bào bất tử hóa phi ác tính. Để thiết lập đối chứng âm thực sự, phải sử dụng dữ liệu phiên mã mô nguyên phát từ **GTEx**.
4. **Tuyệt đối không quy chụp dấu (+/-) từ kết quả Docking:** Docking chỉ trả về Ái lực (Độ lớn). Việc giải mã thuốc kích hoạt hay ức chế BẮT BUỘC phải để Lớp 3 (Dữ liệu RNA-seq / AUCell) quyết định.
5. **Tuyệt đối không bỏ qua "Bể chứa Cạnh tranh":** Phương trình RWR phải sử dụng Trọng số Tích số (Ái lực $\times$ Độ phong phú Protein). Bơm tín hiệu mù quáng mà không biết số lượng bản sao protein trong tế bào sẽ vi phạm Định luật Tác dụng Khối lượng.
6. **Tuyệt đối không để dữ liệu Omics bằng 0 (Omics Dropout):** Khi tính toán $P_0$, nếu dữ liệu biểu hiện gen (TPM) bị khuyết, BẮT BUỘC gán pseudo-count (ví dụ: $E_i = 0.001$). Một số 0 tròn trĩnh sẽ làm sụp đổ phương trình toán học.
7. **Tuyệt đối không để "Nút mồ côi" phá thuật toán:** Trước khi chạy RWR, phải có hàm kiểm tra `Degree > 0`. Nút không có cạnh liên kết phải bị gạt khỏi vector $P_0$ để tránh làm giam cầm xác suất truyền tin.

---

### BẢN THIẾT KẾ WORKFLOW HOÀN CHỈNH (THE MASTER PIPELINE)

#### LAYER 1: CỔNG LỌC NHIỆT ĐỘNG HỌC & MỎ NEO VẬT LÝ (DOCKING)
*Xác định tiềm năng liên kết lý thuyết và thiết lập lực đẩy khởi điểm.*

* **[S] Cụ thể:** Chạy GNINA cho 180 hợp chất trên 9 đích (MMPs, PPARs, COX-2). Chỉ giữ lại các tương tác vượt qua ngưỡng $\Delta G < -6.5$ kcal/mol. Trích xuất chỉ số `CNN_VS`.
* **[M] Đo lường:** Tỷ lệ PoseBusters Pass = 100%. Tạo ra Vector $P_0$ thô.
* **[A] Khả thi:** Tự động hóa bằng bash script và PLIP CLI.
* **[R] Liên quan:** Loại bỏ "Nghịch lý chuẩn hóa $P_0$" ngay từ vòng gửi xe bằng cách chặn đứng các chất có ái lực yếu.
* **[T] Thời hạn:** 1 tuần.
* **FAIR Data:**
    * **[F&A]**: Lưu file `.sdf` theo chuẩn danh pháp `[LigandID]_[TargetPDB].sdf`.
    * **[I&R]**: Ghi log mọi hợp chất bị loại bỏ do không đạt ngưỡng năng lượng vào tệp `L1_Eliminated_Compounds.csv`.

#### LAYER 2: BẢN ĐỒ HÌNH HỌC ĐẶC THÙ (TOPOLOGY & RWR)
*Lan truyền xung lực ức chế dựa trên sự dung hợp giữa không gian vật lý và hệ gen.*

* **[S] Cụ thể:** 1. Lập ma trận bất đối xứng $W_{Integrated}$ bằng cách kết hợp PPI (vô hướng) và pySCENIC (có hướng).
    2. Cập nhật phương trình khởi tạo: $P_{0, i} \propto CNN\_VS_i \times E_i$ (Ái lực $\times$ Biểu hiện gen CCLE) cộng thêm pseudo-count $0.001$.
    3. Chạy RWR bằng `nx.pagerank` với `alpha = 0.7` (tương đương Restart probability $r=0.3$).
* **[M] Đo lường:** Trích xuất chính xác Top 50 Hub Genes chịu xác suất dừng cao nhất.
* **[A] Khả thi:** Khởi chạy tập lệnh Python `layer2_v2.py` với cấu trúc `nx.DiGraph`.
* **[R] Liên quan:** Chuyển đổi Ái lực (Lớp 1) thành Bản đồ Dư chấn (Maximal Topological Footprint), giải quyết triệt để "Nghịch lý Bể chứa Cạnh tranh".
* **[T] Thời hạn:** 2 tuần (tập trung thời gian chạy GRNBoost2).
* **FAIR Data:**
    * **[F&A]**: Xuất mạng lưới dưới dạng `.graphml` để Cytoscape đọc được.
    * **[I&R]**: Ghi rõ version NetworkX, hệ số `alpha` và `seed` vào file `manifest.json`.

#### LAYER 3: MÀNG LỌC DI TRUYỀN & PHIÊN MÃ (VULNERABILITY VALIDATION)
*Cầu dao sinh tử: Giải mã dấu chức năng và thẩm định khả năng thích nghi.*

* **[S] Cụ thể:** 1. Lọc Top 50 Hub Genes qua dữ liệu CRISPR (DepMap). Chỉ giữ lại các gen có Xác suất phụ thuộc (Probability of Dependency) $> 0.8$ tại tế bào ung thư.
    2. Đối chiếu mạng lưới ung thư với mạng lưới khỏe mạnh nguyên phát từ cơ sở dữ liệu **GTEx** (không dùng HEK-293 làm chuẩn khỏe mạnh).
    3. Dùng AUCell của pySCENIC để xác định Hub Gene bị Up-regulated hay Down-regulated.
* **[M] Đo lường:** Số lượng Master Regulons/Hub Genes còn sót lại sau màng lọc. $\Delta$AUC > 0.05 đối với trạng thái phiên mã.
* **[A] Khả thi:** Truy vấn ma trận `CRISPRGeneDependency.csv` bằng Pandas trong vài giây.
* **[R] Liên quan:** * Chống lại "Thiên kiến Xác nhận": Nếu độc tính là do off-target, Lớp 3 sẽ trả về CRISPR $< 0.5$ và hệ thống tự động loại bỏ giả thuyết MoA này.
    * Đập tan "Lệch pha Thời gian": Ngưỡng CRISPR $> 0.8$ chứng minh tế bào ung thư hoàn toàn vô phương cứu chữa, không thể kích hoạt vòng lặp thích nghi trong 72h tới.
* **[T] Thời hạn:** 3 ngày.
* **FAIR Data:**
    * **[F&A]**: Đổi tên cột kết quả thành `Probability_of_Dependency` chuẩn mực.
    * **[I&R]**: Tích hợp danh pháp gen HGNC để máy móc (machine-readable) có thể đối chiếu tự động với GTEx.

#### LAYER 4: HIỆU CHUẨN KIỂU HÌNH (EMPIRICAL DOSE-RESPONSE)
*Sử dụng thực nghiệm in vitro để lấp đầy khoảng trống nồng độ và đột biến cấu trúc.*

* **[S] Cụ thể:** Đặt kết quả đo IC50 (MTT Assay) trên 6 dòng tế bào cạnh dự đoán của Lớp 3.
* **[M] Đo lường:** $IC_{50}$ thực nghiệm và Chỉ số Chọn lọc (Selectivity Index - SI).
* **[A] Khả thi:** Phụ thuộc vào tiến độ phòng lab thực nghiệm.
* **[R] Liên quan:** * Lấp đầy "Khoảng trống Nồng độ": Xác nhận $K_d$ lý thuyết của docking đủ mạnh để vượt qua hàng rào dược động học nội bào.
    * Tự sửa sai Cấu trúc: Nếu tế bào mang đột biến kháng thuốc không được mô phỏng ở Lớp 1, Lớp 4 sẽ trả về $IC_{50}$ cao, tự động đập tan giả thuyết gắn kết.
* **[T] Thời hạn:** Khớp với lịch thí nghiệm in vitro.
* **FAIR Data:**
    * **[F&A]**: Ảnh đồ thị liều-đáp ứng (Dose-response curve) lưu chuẩn `.svg` độ phân giải cao.
    * **[I&R]**: Nén toàn bộ Script, Ma trận kết quả Lớp 1-3, và Raw data Lớp 4 tải lên Zenodo để lấy mã DOI phục vụ công bố Q1.

---
**TỔNG KẾT:** Cỗ máy Dược lý Hệ thống của bạn đã được thiết kế hoàn chỉnh với tính năng tự chẩn đoán, tự sửa lỗi và tích hợp hệ thống kiểm chứng chéo (Popperian Falsification). Hãy bắt tay vào viết dòng code đầu tiên. Không cần lùi bước trước bất kỳ câu hỏi phản biện nào nữa!
