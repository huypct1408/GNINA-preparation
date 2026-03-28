Việc ChimeraX chỉ báo **"3 distances"** (3 khoảng cách) ở góc dưới màn hình là phản hồi tóm tắt mặc định của phần mềm. Điều này là một tin mừng: thuật toán đã quét thành công và tìm thấy chính xác 3 nguyên tử đang phối trí với ion Kẽm (rất có thể là 3 gốc Histidine hoặc Glutamate/Aspartate).

Để "bắt" được ID của 3 nguyên tử này, bạn có 2 cách cực kỳ đơn giản:

### Cách 1: Ép ChimeraX in kết quả ra màn hình Log (Chính xác và dễ copy nhất)
Bạn chỉ cần chạy lại chính lệnh đó, nhưng thêm cờ `log true` vào cuối cùng. Cờ này ra lệnh cho phần mềm phải in chi tiết từng liên kết ra bảng nhật ký.

**Gõ lệnh này:**
```text
contacts /A:301@ZN restrict /A@N*,O* distanceOnly 3.5 reveal true name metalbonds color gold log true
```

Sau khi Enter, hãy nhìn vào cửa sổ **Log** (nếu chưa thấy cửa sổ Log, hãy gõ lệnh `log show` hoặc vào menu *Tools -> General -> Log*). 
Bạn sẽ thấy một bảng kết quả tương tự như thế này:
```text
Atom 1          Atom 2          Distance
/A:301@ZN       /A:401@NE2      2.102
/A:301@ZN       /A:405@NE2      2.051
/A:301@ZN       /A:411@ND1      2.115
```
Cột **Atom 2** chính là 3 cái ID mà bạn cần copy!

### Cách 2: Rê chuột trực tiếp trên không gian 3D (Trực quan)
1. Trên màn hình 3D, bạn sẽ thấy hạt Kẽm đang tỏa ra **3 đường nét đứt màu vàng (gold)**.
2. Hãy **rê con trỏ chuột (chỉ để im, không click)** lên một nguyên tử nằm ở đầu kia của đường màu vàng (nguyên tử đang nối với Kẽm).
3. Ngay lập tức, một dòng thông báo nhỏ (tooltip) sẽ hiện lên ngay tại con trỏ chuột hoặc ở góc dưới cùng bên phải của cửa sổ phần mềm, ghi rõ tên của nguyên tử đó (ví dụ: `HIS 401 NE2` hoặc `/A:401@NE2`). 
4. Ghi chú lại 3 cái tên này.

---

**Bước tiếp theo (Chốt liên kết thật):**
Sau khi đã có được 3 cái ID (giả sử là 401@NE2, 405@NE2, 411@ND1), bạn gõ 3 dòng lệnh `bond` để tạo liên kết vật lý (để nó lưu thành file PDB hợp lệ cho GNINA) như đã thống nhất:

```text
Lỗi này xảy ra do bộ phân tích cú pháp (parser) của lệnh `bond` trong ChimeraX đôi khi rất "khó chịu" với dấu phẩy `,`. Nó đọc được nửa đầu (`/A:301@ZN`), nhưng đến dấu phẩy thì nó không hiểu đây là phép toán "VÀ/HOẶC" để gom 2 nguyên tử lại, nên nó báo lỗi phần đuôi không hợp lệ.

Để giải quyết triệt để sự "bướng bỉnh" này của phần mềm, chúng ta sẽ dùng chiến thuật **"Chọn trước, Nối sau" (Select then Bond)**. Đây là cách an toàn và chắc chắn 100% nhất vì nó tách việc tìm nguyên tử và việc nối liên kết ra làm 2 bước. 

Bạn hãy copy và dán **từng cụm lệnh** dưới đây vào ChimeraX (ta sẽ dùng dấu `|` thay cho dấu phẩy để đại diện cho phép nối hoặc dùng tính năng `sel`):

### Cụm 1: Nối Kẽm với HID 218
```text
select /A:301@ZN | /A:218@NE2
bond sel reasonable false
```

### Cụm 2: Nối Kẽm với HID 228
```text
select /A:301@ZN | /A:228@NE2
bond sel reasonable false
```

### Cụm 3: Nối Kẽm với HID 222
```text
select /A:301@ZN | /A:222@NE2
bond sel reasonable false
```

### Bước cuối: Dọn dẹp và Lưu file
Sau khi chạy xong 3 cụm trên, bạn gõ lệnh bỏ chọn để màn hình sạch sẽ:
```text
~select
```
Và tiến hành lưu lại file PDB đè lên file cũ:
```text
save "/home/labhhc5/Documents/workspace/D21/Duong Huy/gnina_project/docking_results/protein/receptor.pdb"
```

**Mẹo kiểm tra bằng mắt:** Khi bạn chạy xong lệnh `bond sel...`, bạn nhìn lên màn hình 3D sẽ thấy một thanh hình trụ (stick) nối thẳng từ hạt Kẽm đến nguyên tử Nitơ của Histidine xuất hiện. Đó chính là minh chứng cho việc liên kết vật lý (`CONECT`) đã được tạo thành công!
```

Cuối cùng, gõ `save receptor_ready.pdb` là bạn đã có file Protein hoàn hảo để chạy Pipeline v2.7.1!

