// s1102043_FinalProject.cpp : 此檔案包含 'main' 函式。程式會於該處開始執行及結束執行。
//
#include <iostream> 
#include <iomanip>
#include <windows.h>
#include <gdiplus.h>
#include <fstream>
using namespace Gdiplus;
using namespace std;
#pragma comment (lib, "Gdiplus.lib")
void Log(string tag, string s)
{
#ifdef _DEBUG
    cout << "[" << tag << "]" << s << endl;
#endif
}
template<class T>
void fileLog(string filename, T& obj)
{
#ifdef _DEBUG
    ofstream file;
    file.open(filename);
    file << obj;
    file.close();
#endif
}

template <class T>
class MATRIX
{
private:
    int W;
    int H;
    T** A;
public:
    MATRIX() { A = 0; W = H = 0; }
    MATRIX(int h, int w) {
        init(h, w);
    }
    MATRIX(MATRIX& other) {
        init(other.H, other.W);
        *this = other;
    }
    MATRIX(MATRIX&& other) noexcept {
        init(other.H, other.W);
        *this = other;
    }
    ~MATRIX() { clear(); }
    void init(int h, int w)
    {
        W = w;
        H = h;
        A = new T * [H];
        for (int i = 0; i < H; i++)
        {
            A[i] = new T[W];
            for (int j = 0; j < W; j++)
                A[i][j] = 0;
        }
    }
    void clear()
    {
        if (A)
        {
            for (int i = 0; i < H; i++)
                delete[] A[i];
            delete[] A;
        }
        A = 0;
    }
    void set(int i, int j, const T& v)
    {
        A[i][j] = v;
    }

    template <int h, int w>
    void set(T(&m)[h][w]) //將一個二維陣列的值設定到 MATRIX 類別的內部矩陣中
    {
        init(h, w);
        int i;
        int j;
        for (i = 0; i < H; i++)
            for (j = 0; j < W; j++)
                A[i][j] = *(*(m + j) + i);
    }
    void set(MATRIX& other)//將另一個 MATRIX 物件的內容複製到當前的 MATRIX 物件中
    {
        clear();
        init(other.H, other.W);
        int i;
        int j;
        for (i = 0; i < H; i++)
            for (j = 0; j < W; j++)
                A[i][j] = other.A[i][j];
    }
    int getH() { return H; }
    int getW() { return W; }
    int countValue(T v)
    {
        int n = 0;
        for (int y = 0; y < H; y++)
            for (int x = 0; x < W; x++)
                n += (A[y][x] == v);
        return n;
    }

    void bound(T min, T max)//將矩陣中的每個元素限制在指定的範圍內（min 和 max 之間）
    {
        int i;
        int j;
        for (i = 0; i < H; i++) // 遍歷矩陣的每一行
            for (j = 0; j < W; j++) // 遍歷矩陣的每一列
            {
                if (A[i][j] >= max) A[i][j] = max; // 如果元素大於 max，設為 max
                else if (A[i][j] <= min) A[i][j] = min; // 如果元素小於 min，設為 min
            }

    }
    void interbound(T min, T max)//矩陣的範圍縮放操作
    {
        int i;
        int j;
        T a = 0;
        T b = 0;
        a = b = A[0][0];
        for (i = 0; i < H; i++)
            for (j = 0; j < W; j++)
            {
                if (A[i][j] >= a) a = A[i][j];
                else if (A[i][j] <= b) b = A[i][j];
            }
        for (i = 0; i < H; i++)
            for (j = 0; j < W; j++)
            {
                A[i][j] = min + (max - min) * (A[i][j] - b) / (a - b);
            }
    }

    void transpose()//轉置操作
    {
        MATRIX t(W, H);
        for (int i = 0; i < H; i++)
            for (int j = 0; j < W; j++)
                t.A[j][i] = A[i][j];
        clear();
        init(t.H, t.W);
        for (int i = 0; i < H; i++)
            for (int j = 0; j < W; j++)
                A[i][j] = t.A[i][j];

    }
    T dot(MATRIX& m)//內積操作
    {
        T d = 0;
        for (int i = 0; i < H; i++)
            for (int j = 0; j < W; j++)
                d += m.A[i][j] * A[i][j];
        return d;
    }
    void setAll(T c)//將矩陣中的所有元素設置為指定的值 c
    {
        for (int i = 0; i < H; i++)
            for (int j = 0; j < W; j++)
                A[i][j] = c;
    }
    T& operator()(int i, int j) { return A[i - 1][j - 1]; }
    T* operator[](int i) { return A[i]; }

    MATRIX& operator=(const MATRIX& m)
    {
        clear();
        init(m.H, m.W);
        int i;
        int j;
        for (i = 0; i < H; i++)
            for (j = 0; j < W; j++)
                A[i][j] = m.A[i][j];
        return *this;
    }

    void operator*=(const T& s)
    {
        int i;
        int j;
        for (i = 0; i < H; i++)
            for (j = 0; j < W; j++)
                A[i][j] *= s;
    }
    void operator*=(const double& s)
    {
        int i;
        int j;
        for (i = 0; i < H; i++)
            for (j = 0; j < W; j++)
                A[i][j] *= s;
    }
    void operator/=(const T& s)
    {
        int i;
        int j;
        for (i = 0; i < H; i++)
            for (j = 0; j < W; j++)
                A[i][j] /= s;
    }
    void operator+=(const MATRIX<T>& m)
    {
        int i, j;
        for (i = 0; i < H; i++)
            for (j = 0; j < W; j++)
                A[i][j] += m.A[i][j];
    }
    void operator-=(const MATRIX<T>& m)
    {
        int i, j;
        for (i = 0; i < H; i++)
            for (j = 0; j < W; j++)
                A[i][j] -= m.A[i][j];
    }

    MATRIX<T> operator*(const MATRIX<T>& m)
    {
        if (W != m.H) throw exception("can not cross product");
        MATRIX<T> C(H, m.W);
        int i, j, k;
        for (i = 0; i < C.H; i++)
            for (j = 0; j < C.W; j++)
                for (k = 0; k < W; k++)
                    C[i][j] += A[i][k] * m.A[k][j];
        return C;
    }
    MATRIX<T> conv(MATRIX<T> m)//矩陣的卷積操作
    {
        MATRIX<T> C(H, W);
        int i, j;
        int x, y;
        int dx = m.W / 2;
        int dy = m.H / 2;
        int n = m.W * m.H;
        for (y = dy; y < H - dy; y++)
            for (x = dx; x < W - dx; x++)
            {
                C[y][x] = 0;
                for (i = 0; i < m.H; i++)
                    for (j = 0; j < m.W; j++)
                        C[y][x] += (A[y + i - dy][x + j - dx] * m(i + 1, j + 1));
            }
        return C;
    }
    //--------------------------------------------------------
    friend ostream& operator<<(ostream& os, const MATRIX<T>& m)
    {
        int i;
        int j;
        for (i = 0; i < m.H; i++)
        {
            for (j = 0; j < m.W; j++)
                os << setw(5) << m.A[i][j];//setw(5)為將下一個輸出的值寬度設為 5
            os << endl;
        }
        return os;
    }
    friend istream& operator>>(istream& is, MATRIX<T>& m)
    {
        int i;
        int j;
        for (i = 0; i < m.H; i++)
            for (j = 0; j < m.W; j++)
                is >> m.A[i][j];
        return is;
    }
};

class APP
{
private:
    Gdiplus::GdiplusStartupInput gdiplusStartupInput;
    ULONG_PTR gdiplusToken;
    Gdiplus::Status status;
    MATRIX<int> A[3];
    int W = 0;
    int H = 0;
    std::wstring wfilename;
    Bitmap* bmp = 0;
	
    int GetEncoderClsid(const WCHAR* format, CLSID* pClsid) // 獲取指定格式的編碼器 CLSID
    {
        UINT num = 0, size = 0;
        GetImageEncodersSize(&num, &size);
        if (size == 0) return -1;

        ImageCodecInfo* pImageCodecInfo = (ImageCodecInfo*)(malloc(size));
        if (!pImageCodecInfo) return -1;

        GetImageEncoders(num, size, pImageCodecInfo);
        for (UINT j = 0; j < num; ++j) {
            if (wcscmp(pImageCodecInfo[j].MimeType, format) == 0) {
                *pClsid = pImageCodecInfo[j].Clsid;
                free(pImageCodecInfo);
                return j;
            }
        }
        free(pImageCodecInfo);
        return -1;
    }

public:
    APP()
    {
		status = Gdiplus::GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);// 初始化 GDI+ 庫
    }
    ~APP()
    {
        if (bmp) delete bmp;
        Gdiplus::GdiplusShutdown(gdiplusToken);
    }
    int getH() { return H; }
    int getW() { return W; }

    bool read(string filename)
    {
        Log("bitmap", "read");
        wfilename = std::wstring(filename.begin(), filename.end());
        bmp = Bitmap::FromFile(wfilename.c_str());
        W = bmp->GetWidth();
        H = bmp->GetHeight();
        Color pixelColor;
        for (int i = 0; i < 3; i++)
            A[i].init(H, W);
        for (int y = 0; y < H; y++) {
            for (int x = 0; x < W; x++) {
                bmp->GetPixel(x, y, &pixelColor);
                A[2](y + 1, x + 1) = pixelColor.GetRed();
                A[1](y + 1, x + 1) = pixelColor.GetGreen();
                A[0](y + 1, x + 1) = pixelColor.GetBlue();
            }
        }
        return true;
    }
    bool write(string filename)
    {
        Log("bitmap", "write");
        wstring wfilename = std::wstring(filename.begin(), filename.end());
        CreateAndSaveBitmap(wfilename, PixelFormat24bppRGB, W, H, L"image/bmp");
        return true;
    }
    CLSID GetEncoderClsid(const WCHAR* format) {
        UINT num = 0, size = 0;
		GetImageEncodersSize(&num, &size);//獲取編碼器數量和大小
        if (size == 0) return CLSID{};
        ImageCodecInfo* pImageCodecInfo = (ImageCodecInfo*)(malloc(size));
        if (pImageCodecInfo == nullptr) return CLSID{};
        GetImageEncoders(num, size, pImageCodecInfo);
        for (UINT i = 0; i < num; ++i) {
            if (wcscmp(pImageCodecInfo[i].MimeType, format) == 0) {
                CLSID clsid = pImageCodecInfo[i].Clsid;
                free(pImageCodecInfo);
                return clsid;
            }
        }
        free(pImageCodecInfo);
        return CLSID{};
    }

    void CreateAndSaveBitmap(const std::wstring& filename, PixelFormat format, int width, int height, const std::wstring& mimeType)
    {
        Bitmap bmp(width, height, format);
        Rect rect(0, 0, width, height);
        BitmapData bmpData;
        bmp.LockBits(&rect, ImageLockModeWrite, format, &bmpData);
        BYTE* scan0 = static_cast<BYTE*>(bmpData.Scan0);
        int stride = bmpData.Stride;
        int bytesPerPixel = (format == PixelFormat24bppRGB) ? 3 : 4;
        for (int y = 0; y < height; ++y) {
            BYTE* row = scan0 + y * stride;
            for (int x = 0; x < width; ++x) {
                BYTE* pixel = row + x * bytesPerPixel;
                pixel[0] = A[0](y + 1, x + 1);
                pixel[1] = A[1](y + 1, x + 1);
                pixel[2] = A[2](y + 1, x + 1);
            }
        }
        bmp.UnlockBits(&bmpData);
        CLSID clsid = GetEncoderClsid(mimeType.c_str());
        if (clsid == CLSID{}) {
            std::wcerr << L"Failed to get encoder CLSID for " << mimeType << std::endl;
            return;
        }
        Status stat = bmp.Save(filename.c_str(), &clsid, nullptr);
        if (stat == Ok)
            std::wcout << L"Saved: " << filename << std::endl;
        else
            std::wcerr << L"Save failed: " << stat << std::endl;
    }
    void updateImage()
    {
        Color pixelColor;

        for (int y = 0; y < H; y++) {
            for (int x = 0; x < W; x++) {
                pixelColor.SetValue(pixelColor.MakeARGB(255, A[0](y + 1, x + 1), A[1](y + 1, x + 1), A[2](y + 1, x + 1)));
                bmp->SetPixel(x, y, pixelColor);
            }
        }
    }
    void toBiLevel()//二值化（灰階轉黑白）
    {
        Log("bitmap", "to bilevel");
        int y;
        int x;
        for (y = 0; y < H; y++)
            for (x = 0; x < W; x++)
            {
                int c = int(0.299 * A[0](y + 1, x + 1) + 0.587 * A[1](y + 1, x + 1) + 0.114 * A[2](y + 1, x + 1));

                if (c > 128) A[0](y + 1, x + 1) = A[1](y + 1, x + 1) = A[2](y + 1, x + 1) = 255;
                else  A[0](y + 1, x + 1) = A[1](y + 1, x + 1) = A[2](y + 1, x + 1) = 0;
            }

    }
    void diffX(int i = 0)//邊緣檢測（亮度差異）
    {
        Log("bitmap", "diffX");
        int y;
        int x;
        for (y = 0; y < H; y++)
            for (x = 0; x < W - 1; x++)
                if ((A[i][y][x] - A[i][y][x + 1]) * (A[i][y][x] - A[i][y][x + 1]) > 50) A[i][y][x] = 0;
                else A[i][y][x] = 255;
    }
    void diffY(int i = 0)//邊緣檢測（亮度差異）
    {
        Log("bitmap", "diffY");
        int y;
        int x;
        for (y = 0; y < H - 1; y++)
            for (x = 0; x < W; x++)
                if ((A[i][y][x] - A[i][y + 1][x]) * (A[i][y][x] - A[i][y + 1][x]) > 50) A[i][y][x] = 0;
                else A[i][y][x] = 255;
    }
    void diff(int i = 0)//邊緣檢測（亮度差異）
    {
        Log("bitmap", "diff");
        int y;
        int x;
        for (y = 0; y < H - 1; y++)
            for (x = 0; x < W - 1; x++)
                if ((A[i][y][x] - A[i][y][x + 1]) * (A[i][y][x] - A[i][y][x + 1]) > 50 || (A[i][y][x] - A[i][y + 1][x]) * (A[i][y][x] - A[i][y + 1][x])) A[i][y][x] = 0;
                else A[i][y][x] = 255;
    }
    void adjLight(int i, int delta)//調整亮度（加亮或變暗）
    {
        Log("bitmap", "darker");
        MATRIX<int> diff;
        diff = A[i];
        diff.setAll(delta);
        A[i] += diff;
    }
    void neg()//負片效果
    {
        Log("bitmap", "neg");
        MATRIX<int> m[3];
        for (int i = 0; i < 3; i++)
        {
            m[i] = A[i];
            m[i].setAll(255);
            m[i] -= A[i];
            A[i] = m[i];
        }
    }
    void sampling(double r)//降採樣處理
    {
        Log("bitmap", "sampling");
        int y;
        int x;
        MATRIX<int> B[3];// 用來存放降採樣後的影像數據
        for (int i = 0; i < 3; i++)
        {
            B[i] = A[i];
            B[i].setAll(255);
            //-------------------------------------------
            for (y = 0; y < H / r; y++)
                for (x = 0; x < W / r; x++)
                    B[i][y][x] = A[i][int(y * r)][int(x * r)];// 根據比例 r，從原始影像中選取對應的像素值，填入降採樣後的矩陣
            //-------------------------------------------
            A[i] = B[i];
        }
    }
    void mirror()//鏡射
    {
        Log("bitmap", "mirror");
        int y;
        int x;
        MATRIX<int> B[3];
        for (int i = 0; i < 3; i++)
        {
            B[i].init(H, W);
            //-------------------------------------------
            for (y = 0; y < H; y++)
                for (x = 0; x < W; x++)
                    B[i](y + 1, x + 1) = A[i](H - y, x + 1);
            //-------------------------------------------
            A[i] = B[i];
        }
    }
    void mirror_diagonal()//對角線鏡射操作
    {
        Log("bitmap", "diagonal mirror");
        for (int i = 0; i < 3; i++)
        {
            A[i].transpose();
            W = A[i].getW();
            H = A[i].getH();
        }
        bmp->RotateFlip(Rotate90FlipNone);
    }
	void rotate(int deg = 0)//旋轉
    {
        Log("bitmap", "rotate");
        int y;
        int x;
        MATRIX<int> B[3];
        for (int i = 0; i < 3; i++)
        {
            B[i] = A[i];
            B[i].setAll(255);
            //-------------------------------------------
            for (y = 0; y < H; y++)
                for (x = 0; x < W; x++)
                {
                    double theta = deg * 3.14159 / 180.;
                    int a = x * cos(theta) + y * sin(theta);
                    int b = x * sin(theta) - y * cos(theta);
                    if (a >= 0 && a < H && b >= 0 && b < W)
                        B[i][y][x] = A[i][a][b];
                }
            //-------------------------------------------
            A[i] = B[i];
        }
    }
	void offset(int dx = 0, int dy = 0)//位移操作
    {
        Log("bitmap", "offset");
        int y;
        int x;
        MATRIX<int> B[3];
        for (int i = 0; i < 3; i++)
        {
            B[i] = A[i];
            B[i].setAll(255);
            //-------------------------------------------
            for (y = 0; y < H; y++)
                for (x = 0; x < W; x++)
                {
                    int a = y + dy;
                    int b = x + dx;
                    if (a >= 0 && a < H && b >= 0 && b < W)
                        B[i][y][x] = A[i][a][b];
                }
            //-------------------------------------------
            A[i] = B[i];
        }
    }
	void avg()//平均化處理
    {
        Log("bitmap", "avg");
        int m[3][3] = { {1,1,1},{1,1,1},{1,1,1} };
        MATRIX<int> mask;
        mask.set(m);
        MATRIX<int> B[3];
        for (int i = 0; i < 3; i++)
        {
            B[i] = A[i];
            A[i] = B[i].conv(mask);
            A[i] /= 9;
        }
    }
	void sobel()// Sobel 邊緣檢測
    {
        Log("bitmap", "sobel");
        int m[3][3] = { {-1,0,1},{-2,0,2},{-1,0,1} };
        MATRIX<int> mask;
        mask.set(m);
        MATRIX<int> B[3];
        for (int i = 0; i < 3; i++)
        {
            B[i] = A[i];
            A[i] = B[i].conv(mask);
            A[i].interbound(0, 255);
        }

    }
	void laplacian()// 拉普拉斯邊緣檢測
    {
        Log("bitmap", "laplacian");
        int m[3][3] = { {0, -1, 0} ,{ -1, 5, -1},{0, -1, 0} };
        MATRIX<int> mask;
        mask.set(m);
        MATRIX<int> B[3];
        for (int i = 0; i < 3; i++)
        {
            B[i] = A[i];
            A[i] = B[i].conv(mask);
            A[i].interbound(0, 255);
        }

    }

    void histogramEqualization()
    {
        for (int i = 0; i < 3; i++) // 對 RGB 三個色彩通道分別進行處理
        {
            int histogram[256] = { 0 }; // 儲存每個像素值（0-255）的出現次數
            int cdf[256] = { 0 };       // 累積分佈函數（Cumulative Distribution Function）
            int mapping[256] = { 0 };   // 像素值的映射表
            int totalPixels = H * W;    // 總像素數量

            // 1. 計算直方圖
            for (int y = 0; y < H; y++)
                for (int x = 0; x < W; x++)
                    histogram[A[i][y][x]]++;

            // 2. 計算累積分佈函數（CDF）
            cdf[0] = histogram[0];
            for (int j = 1; j < 256; j++)
                cdf[j] = cdf[j - 1] + histogram[j];

            // 3. 計算像素值的映射表
            for (int j = 0; j < 256; j++)
                mapping[j] = int((float)cdf[j] / totalPixels * 255);

            // 4. 將原始像素值映射到新的像素值
            for (int y = 0; y < H; y++)
                for (int x = 0; x < W; x++)
                    A[i][y][x] = mapping[A[i][y][x]];
        }

        updateImage(); // 更新影像數據
        write("equalized.bmp"); // 將均衡化後的影像儲存為檔案
    }



    void show()
    {
        HWND hWnd = CreateWindow(L"edit", L"wp", WS_OVERLAPPEDWINDOW, 0, 0, bmp->GetWidth(), bmp->GetHeight(), GetDesktopWindow(), (HMENU)0, 0, 0);
        ShowWindow(hWnd, SW_SHOW);
        HDC hdc = GetDC(hWnd);
        Graphics graphics(hdc);
        CachedBitmap cachedBitmap(bmp, &graphics);
        graphics.DrawCachedBitmap(&cachedBitmap, 0, 0);
        cin.ignore();
        DestroyWindow(hWnd);
    }
    //-------------------------------------------
	void pixelProcess()//像素處理
    {
        Log("pixelProcess", "------------------------");
        read("logo.bmp");
        write("logo_gray.bmp");
        adjLight(0, -20);
        write("logo_darker.bmp");
        toBiLevel();
        write("logo_bilevel.bmp");
        neg();
        write("logo_neg.bmp");
        diff();
        write("logo_diff.bmp");
    }
	void geometryProcess()//幾何處理
    {
        Log(" geometryProcess", "------------------------");
        read("logo.bmp");
        write("logo_gray.bmp");
        mirror();
        write("logo_mirror.bmp");
        offset(300, 300);
        write("logo_offset.bmp");
        sampling(2);
        write("logo_sampling.bmp");
        rotate(60);
        write("logo_rotate.bmp");
        mirror_diagonal();
        write("logo_mirror_diagonal.bmp");
    }
	void convProcess()//卷積處理
    {
        Log(" convProcess", "------------------------");
        read("logo.bmp");
        write("logo_gray.bmp");
        avg();
        write("logo_avg.bmp");
        laplacian();
        write("logo_laplacian.bmp");
        sobel();
        write("logo_sobel.bmp");
    }
	void fusion()//融合處理
    {
        Log(" hybridProcess", "fusion");
        read("logo.bmp");
        int m[3][3] = { {-1,0,1},{-2,0,2},{-1,0,1} };
        MATRIX<int> mask;
        mask.set(m);
        MATRIX<int> B[3];
        for (int i = 0; i < 3; i++)
        {
            B[i] = A[i];
            A[i] += B[i].conv(mask);
            A[i].bound(0, 255);
        }
        updateImage();
        write("logo_fusion.bmp");
        // show();
    }
	void substract()//減法處理
    {
        Log(" hybridProcess", "substract");
        read("logo.bmp");
        MATRIX<int> B[3];
        for (int i = 0; i < 3; i++)
        {
            B[i] = A[i];
            offset(1, 1);
            A[i] -= B[i];
            A[i].bound(0, 255);
            neg();
            A[i].bound(0, 255);
        }
        updateImage();
        write("logo_substract.bmp");
    }
	void hybridProcess()//混合處理
    {
        Log(" hybridProcess", "------------------------");
        fusion();
        substract();
    }
	void bitplane(string filename, int channel, int i)//位平面處理
    {
        Log("bitmap", "bitplaned");
        read(filename);
        int y;
        int x;
        for (y = 0; y < H; y++)
            for (x = 0; x < W; x++)
            {
                // 將其他兩個色彩通道的值設為 0
                A[(channel + 1) % 3][y][x] = 0;
                A[(channel + 2) % 3][y][x] = 0;

                // 提取指定色彩通道的第 i 位平面
                A[channel][y][x] = ((A[channel][y][x] >> i) & 0x1) * 255;

                // 將提取的位平面值複製到其他兩個色彩通道，生成灰階影像
                A[(channel + 1) % 3][y][x] = A[(channel + 2) % 3][y][x] = A[channel][y][x];
            }

    }
	void bitplanes(string filename, int channel = 0)//位平面處理
    {
        Log(" bitplane", "------------------------");
        int i;
        for (i = 0; i < 8; i++)// 遍歷 8 個位平面
        {
            bitplane(filename.c_str(), channel, i); // 提取第 i 個位平面
            char f[256] = { 0 };
            sprintf_s(f, "bitplane%d.bmp", i);// 生成對應的檔案名稱，例如 "bitplane0.
            write(f);
        }
    }
	void datahiding_encode()//數據隱藏編碼
    {
        Log(" data hiding", "encode");
        int y;
        int x;

        read("secret.bmp");
        toBiLevel();

        MATRIX<int> B[3];
        for (int i = 0; i < 3; i++)
        {
            B[i] = A[i];
            for (y = 0; y < H; y++)
                for (x = 0; x < W; x++)
                    B[i][y][x] = ((B[i][y][x] & 0xff) >> 5);
            read("logo.bmp");
            for (y = 0; y < H; y++)
                for (x = 0; x < W; x++)
                    A[i][y][x] = (A[i][y][x] & 0xFB);

            for (y = 0; y < H; y++)
                for (x = 0; x < W; x++)
                    A[i][y][x] = A[i][y][x] | B[i][y][x];
            write("stegoimage.bmp");
        }
    } 
	void datahiding_decode()//數據隱藏解碼
    {
        Log(" data hiding", "decode");
        bitplane("stegoimage.bmp", 0, 2);
        write("recoverimage.bmp");
    }
	void datahiding()//數據隱藏處理
    {
        Log(" data hiding", "------------------------");
        datahiding_encode();
        datahiding_decode();

    }

    //--------------------Filter濾鏡------------------------//

    void Comic_Style_Filter()//漫畫風濾鏡
    {
        Log("Fliter", "Use  ComicStyle Fliter!!");
        MATRIX<int> eage[3];
        for (int i = 0; i < 3; i++)
        {
            eage[i] = A[i];
        }
        //邊緣偵測
        int sobelKernel[3][3] = { {-1,0,1},{-2,0,2},{-1,0,1} };
        MATRIX<int> mask;
        mask.set(sobelKernel);
        for (int i = 0; i < 3; i++)
        {
            eage[i] = A[i].conv(mask);
            eage[i].interbound(0, 255);
        }
        //顏色量化
        for (int i = 0; i < 3; i++)
        {
            for (int y = 0; y < H; y++)
            {
                for (int x = 0; x < W; x++)
                {
                    int v = A[i][y][x];
                    A[i][y][x] = (v / 64) * 64;
                }
            }
        }

        //使用eage讓邊緣變黑線
        for (int y = 0; y < H; y++)
        {
            for (int x = 0; x < W; x++)
            {
                int intentstify = (eage[0][y][x] + eage[1][y][x] + eage[2][y][x]) / 3;
                if (intentstify < 100)
                {
                    eage[0][y][x] = eage[1][y][x] = eage[2][y][x] = 0;
                }
            }
        }

        updateImage();
        write("Comic_Style_Image.bmp");

    }

    void OilPainting_Style_Filter()//油畫風濾鏡
    {

    }

    void Mosaic_Style_Filter(int blocksize=10)//馬賽克濾鏡
    {
        Log("Fliter", "Use MosaicStyle Fliter!!");
        for (int i = 0; i < 3; i++)
        {
            for (int y = 0; y < H; y += blocksize)
            {
                for (int x = 0; x < W; x += blocksize)
                {
                    int sum = 0;
                    int count = 0;
                    for (int dy = 0; dy < blocksize; dy++)
                    {
                        for (int dx = 0; dx < blocksize; dx++)
                        {
                            int xx = x + dx;
                            int yy = y + dy;

                            if (yy < H && xx < W)
                            {
                                sum += A[i][yy][xx];
                                count++;
                            }
                        }
                    }
                    int avg = count > 0 ? sum / count : 0;
                    for (int dy = 0; dy < blocksize; dy++)
                    {
                        for (int dx = 0; dx < blocksize; dx++)
                        {
                            int xx = x + dx;
                            int yy = y + dy;

                            if (yy < H && xx < W)
                            {
                                A[i][yy][xx] = avg;
                            }
                           
                        }
                    }
                }
            }
        }
        updateImage();
        write("Mosaic_Style_Image.bmp");
    }
    bool run()
    {
        /*if (status == Gdiplus::Ok) {
             pixelProcess();
             geometryProcess();
             convProcess();
             bitplanes("logo.bmp");
             datahiding();
             return true;
         }
         else
             return false;*/

        if (status == Gdiplus::Ok) {
            read("logo.bmp");
            Comic_Style_Filter();            
            return true;
        }
        else {
            return false;
        }
    }

};

int main() {
    APP app;
    app.run();
	return 0;//0代表程式成功結束，1代表程式異常結束
}