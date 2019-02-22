#include"gdal_priv.h"
#include"cpl_conv.h" // for CPLMalloc()
#include<iostream>
#include<fstream>
#include<string>
#include"Eigen/Dense"

using namespace std;
using namespace Eigen;

double max(double &a, double &b)
{
	if (a > b)
		return a;
	else
		return b;
}

double min(double &a, double &b)
{
	if (a < b)
		return a;
	else
		return b;
}

int main()
{

	GDALAllRegister();          //GDAL所有操作都需要先注册格式
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");  //支持中文路径
	const char* imgPath1 = "D:/gdalData/match_mosaic/50051016_0_rec.tif";
	const char* imgPath2 = "D:/gdalData/match_mosaic/50051017_0_rec.tif";

	GDALDataset* poDataset1 = (GDALDataset *)GDALOpen(imgPath1, GA_ReadOnly);
	if (poDataset1 == nullptr)
	{
		cout << "Can't Open Image!" << endl;
		return 1;
	}
	GDALRasterBand *poBand11, *poBand12, *poBand13;
	poBand11 = poDataset1->GetRasterBand(1);
	poBand12 = poDataset1->GetRasterBand(2);
	poBand13 = poDataset1->GetRasterBand(3);
	
	int imgWidth1 = poDataset1->GetRasterXSize();   //图像宽度
	int imgHeight1 = poDataset1->GetRasterYSize();  //图像高度
	int bandNum1 = poDataset1->GetRasterCount();    //波段数
	int depth1 = GDALGetDataTypeSize(poDataset1->GetRasterBand(1)->GetRasterDataType())/8;    //图像深度
	double adfGeoTransform1[6];
	poDataset1->GetGeoTransform(adfGeoTransform1);
	double minX1 = adfGeoTransform1[0];
	double maxY1 = adfGeoTransform1[3];
	double pixelWidth1 = adfGeoTransform1[1];
	double pixelHeight1 = adfGeoTransform1[5];
	double maxX1 = minX1 + (imgWidth1 * pixelWidth1);
	double minY1 = maxY1 + (imgHeight1 * pixelHeight1);

	GDALDataset* poDataset2 = (GDALDataset *)GDALOpen(imgPath2, GA_ReadOnly);
	if (poDataset2 == nullptr)
	{
		cout << "Can't Open Image!" << endl;
		return 1;
	}
	GDALRasterBand *poBand21, *poBand22, *poBand23;
	poBand21 = poDataset2->GetRasterBand(1);
	poBand22 = poDataset2->GetRasterBand(2);
	poBand23 = poDataset2->GetRasterBand(3);
	
	int imgWidth2 = poDataset2->GetRasterXSize();   //图像宽度
	int imgHeight2 = poDataset2->GetRasterYSize();  //图像高度
	int bandNum2 = poDataset2->GetRasterCount();    //波段数
	int depth2 = GDALGetDataTypeSize(poDataset2->GetRasterBand(1)->GetRasterDataType())/8;    //图像深度
	double adfGeoTransform2[6];
	poDataset2->GetGeoTransform(adfGeoTransform2);
	double minX2 = adfGeoTransform2[0];
	double maxY2 = adfGeoTransform2[3];
	double pixelWidth2 = adfGeoTransform2[1];
	double pixelHeight2 = adfGeoTransform2[5];
	double maxX2 = minX2 + (imgWidth2 * pixelWidth2);
	double minY2 = maxY2 + (imgHeight2 * pixelHeight2);


	double minX = min(minX1, minX2);
	double maxX = max(maxX1, maxX2);
	double minY = min(minY1, minY2);
	double maxY = max(maxY1, maxY2);

	int cols = int((maxX - minX) / pixelWidth1);
	int rows = int((maxY - minY) / abs(pixelHeight1));

	int xOffset1 = int((minX1 - minX) / pixelWidth1);
	int yOffset1 = int((maxY1 - maxY) / pixelHeight1);

	int xOffset2 = int((minX2 - minX) / pixelWidth1);
	int yOffset2 = int((maxY2 - maxY) / pixelHeight1);

	GDALDriver *pDriver = GetGDALDriverManager()->GetDriverByName("GTIFF"); //图像驱动
	char** ppszOptions = NULL;
	
	const char* dstPath = "D:/tts3ixi3.tif";
	
	GDALDataset* dst = pDriver->Create(dstPath, cols, rows, 3, GDT_Byte, ppszOptions);
	if (dst == nullptr)
	{
		printf("Can't Write Image!");
		return false;
	}
	//double adfGeoTransform1[6];
	adfGeoTransform1[0] = adfGeoTransform2[0];
	dst->SetGeoTransform(adfGeoTransform1);//设置坐标

	//申请buf
	size_t imgBufNum = (size_t)cols * rows * depth1;
	GByte *imgBuf1 = new GByte[imgBufNum];
	GByte *imgBuf2 = new GByte[imgBufNum];
	GByte *imgBuf3 = new GByte[imgBufNum];

	size_t imgBufNum1 = (size_t)imgWidth1 * imgHeight1 * 1 * depth1;
	size_t imgBufNum2 = (size_t)imgWidth2 * imgHeight2 * 1 * depth1;
	GByte *imgBuf11 = new GByte[imgBufNum1];
	GByte *imgBuf12 = new GByte[imgBufNum1];
	GByte *imgBuf13 = new GByte[imgBufNum1];
	GByte *imgBuf21 = new GByte[imgBufNum2];
	GByte *imgBuf22 = new GByte[imgBufNum2];
	GByte *imgBuf23 = new GByte[imgBufNum2];


	poBand11->RasterIO(GF_Read, 0, 0, imgWidth1, imgHeight1, imgBuf11, imgWidth1, imgHeight1,
		GDT_Byte, 0, 0);
	poBand12->RasterIO(GF_Read, 0, 0, imgWidth1, imgHeight1, imgBuf12, imgWidth1, imgHeight1,
		GDT_Byte, 0, 0);
	poBand13->RasterIO(GF_Read, 0, 0, imgWidth1, imgHeight1, imgBuf13, imgWidth1, imgHeight1,
		GDT_Byte, 0, 0);

	for (int i = 0; i < imgBufNum; i++)
		imgBuf1[i] = 0;
	for (int i = 0; i < imgBufNum; i++)
		imgBuf2[i] = 0;
	for (int i = 0; i < imgBufNum; i++)
		imgBuf3[i] = 0;
	for (int i = 0; i < yOffset2; i++)
	{
		for (int j = xOffset1; j < cols; j++)
			imgBuf1[i*cols + j] = imgBuf11[i*imgWidth1 + (j - xOffset1)];
	}
	for (int i = 0; i < yOffset2; i++)
	{
		for (int j = xOffset1; j < cols; j++)
			imgBuf2[i*cols + j] = imgBuf12[i*imgWidth1 + (j - xOffset1)];
	}
	for (int i = 0; i < yOffset2; i++)
	{
		for (int j = xOffset1; j < cols; j++)
			imgBuf3[i*cols + j] = imgBuf13[i*imgWidth1 + (j - xOffset1)];
	}


	
	poBand21->RasterIO(GF_Read, 0, 0, imgWidth2, imgHeight2, imgBuf21, imgWidth2, imgHeight2,
		GDT_Byte, 0, 0);
	poBand22->RasterIO(GF_Read, 0, 0, imgWidth2, imgHeight2, imgBuf22, imgWidth2, imgHeight2,
		GDT_Byte, 0, 0);
	poBand23->RasterIO(GF_Read, 0, 0, imgWidth2, imgHeight2, imgBuf23, imgWidth2, imgHeight2,
		GDT_Byte, 0, 0);


	//重叠区域读取并平滑

	for (int i = yOffset2; i < imgHeight1 ; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (j <= cols - imgWidth1)
				imgBuf1[i*cols + j] = imgBuf21[(i - yOffset2)*imgWidth2 + j];
			if (j >(cols - imgWidth1) && j <= imgWidth2)
			{
				if (imgBuf21[(i - yOffset2)*imgWidth2 + j] == 0 || imgBuf11[i *imgWidth1 + j - (cols - imgWidth1)] == 0)
					imgBuf1[i*cols + j] = imgBuf21[(i - yOffset2)*imgWidth2 + j] + imgBuf12[i *imgWidth1 + j - (cols - imgWidth1)];
				else
				{
					imgBuf1[i*cols + j] = (imgBuf21[(i - yOffset2)*imgWidth2 + j] + imgBuf11[i *imgWidth1 + j - (cols - imgWidth1)]) / 2;

				}
				//imgBuf1[i*cols + j] = (imgBuf21[(i - yOffset2)*imgWidth2 + j] + imgBuf11[i *imgWidth1 + j - (cols - imgWidth1)]) / 2;
				//imgBuf1[i*cols + j] = 0;
			}
			if (j > imgWidth2&&j < cols)
				imgBuf1[i*cols + j] = imgBuf11[i*imgWidth1 + j-(cols - imgWidth1)];
		}
	}

	for (int i = yOffset2; i < imgHeight1; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (j <= cols - imgWidth1)
				imgBuf2[i*cols + j] = imgBuf22[(i - yOffset2)*imgWidth2 + j];
			if (j >(cols - imgWidth1) && j <= imgWidth2)
			{ 
				if(imgBuf22[(i - yOffset2)*imgWidth2 + j]==0 || imgBuf12[i *imgWidth1 + j - (cols - imgWidth1)] == 0)
					imgBuf2[i*cols + j] = imgBuf22[(i - yOffset2)*imgWidth2 + j] + imgBuf12[i *imgWidth1 + j - (cols - imgWidth1)];
				else
				{
					imgBuf2[i*cols + j] = (imgBuf22[(i - yOffset2)*imgWidth2 + j] + imgBuf12[i *imgWidth1 + j - (cols - imgWidth1)]) / 2;

				}
				//imgBuf2[i*cols + j] = (imgBuf22[(i - yOffset2)*imgWidth2 + j] + imgBuf12[i *imgWidth1 + j-(cols - imgWidth1)])/2;
				//imgBuf2[i*cols + j] = 0;
			}
			if (j > imgWidth2&&j < cols)
				imgBuf2[i*cols + j] = imgBuf12[i*imgWidth1 + j-(cols - imgWidth1)];
		}
	}

	for (int i = yOffset2; i < imgHeight1; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (j <= cols - imgWidth1)
				imgBuf3[i*cols + j] = imgBuf23[(i - yOffset2)*imgWidth2 + j];
			if (j >(cols - imgWidth1) && j <= imgWidth2)
			{
				if (imgBuf23[(i - yOffset2)*imgWidth2 + j] == 0 || imgBuf13[i *imgWidth1 + j - (cols - imgWidth1)] == 0)
					imgBuf3[i*cols + j] = imgBuf23[(i - yOffset2)*imgWidth2 + j] + imgBuf13[i *imgWidth1 + j - (cols - imgWidth1)];
				else
				{
					imgBuf3[i*cols + j] = (imgBuf23[(i - yOffset2)*imgWidth2 + j] + imgBuf13[i *imgWidth1 + j - (cols - imgWidth1)]) / 2;

				}
				//imgBuf2[i*cols + j] = (imgBuf22[(i - yOffset2)*imgWidth2 + j] + imgBuf12[i *imgWidth1 + j-(cols - imgWidth1)])/2;
				//imgBuf2[i*cols + j] = 0;
			}
			if (j > imgWidth2&&j < cols)
				imgBuf3[i*cols + j] = imgBuf13[i*imgWidth1 + j - (cols - imgWidth1)];
		}
	}
















	for (int i = imgHeight1; i < rows; i++)
	{
		for (int j = 0; j < imgWidth2; j++)
			imgBuf1[i*cols + j] = imgBuf21[(i - imgHeight1+imgHeight1- yOffset2)*imgWidth2 + j];
	}
	 
	for (int i = imgHeight1; i < rows; i++)
	{
		for (int j = 0; j < imgWidth2; j++)
			imgBuf2[i*cols + j] = imgBuf22[(i - imgHeight1 + imgHeight1 - yOffset2)*imgWidth2 + j];
	}
	for (int i = imgHeight1; i < rows; i++)
	{
		for (int j = 0; j < imgWidth2; j++)
			imgBuf3[i*cols + j] = imgBuf23[(i - imgHeight1 + imgHeight1 - yOffset2)*imgWidth2 + j];
	}

	/*for (int i = yOffset2; i < rows; i++)
	{
		for (int j = 0; j < imgWidth2; j++)
			imgBuf1[i*cols + j] = imgBuf21[(i- yOffset2)*imgWidth2 + j];
	}
	for (int i = yOffset2; i < rows; i++)
	{
		for (int j = 0; j < imgWidth2; j++)
			imgBuf2[i*cols + j] = imgBuf22[(i - yOffset2)*imgWidth2 + j];
	}
	for (int i = yOffset2; i < rows; i++)
	{
		for (int j = 0; j < imgWidth2; j++)
			imgBuf3[i*cols + j] = imgBuf23[(i - yOffset2)*imgWidth2 + j];
	}*/

	//delete[] imgBuf21;


	//写入
	
	dst->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, cols, rows, imgBuf1, cols, rows, GDT_Byte, 0, 0);
	dst->GetRasterBand(2)->RasterIO(GF_Write, 0, 0, cols, rows, imgBuf2, cols, rows, GDT_Byte, 0, 0);
	dst->GetRasterBand(3)->RasterIO(GF_Write, 0, 0, cols, rows, imgBuf3, cols, rows, GDT_Byte, 0, 0);

	/*dst->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, cols, rows,
		imgBuf, cols , rows , GDT_Byte, 0, 0);*/
	////释放
	/*delete[] imgBuf;
	imgBuf = nullptr;*/
	//GDALClose(dst);
	
	
	//return 0;

}