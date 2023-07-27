Data = readmatrix("Exper_3_file.xlsx", "Range", "ET4:GD147");

% 데이터 가공
dates = Data(:, 1);
sectorIndices = [2 3; 4 5; 6 7; 8 9; 10 11; 12 13; 14 15; 16 17; 18 19; 20 21; ...
                 22 23; 24 25; 26 27; 28 29; 30 31; 32 33; 34 35; 36 37];

% 폴더 생성
folderName = 'data';
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

% CSV 파일로 저장
fileNames = ["Total_data.csv", "Mine_data.csv", "Manufacturing_data.csv", "Electric_data.csv", ...
             "Waterworks_data.csv", "Construction_data.csv", "Wholesale_data.csv", ...
             "Transportation_data.csv", "Accomodation_data.csv", "IT_data.csv", ...
             "Finance_data.csv", "Real Estate_data.csv", "Scientific Service_data.csv", ...
             "Business Facility management_data.csv", "Education Service_data.csv", ...
             "Health and Social welfare services_data.csv", "Arts Sports Service_data.csv", ...
             "Repair and Other personal services_data.csv"];

for i = 1:numel(fileNames)
    sectorData = [dates, Data(:, sectorIndices(i,2)), Data(:, sectorIndices(i,1))];
    csvwrite(fileNames(i), sectorData);
end

%% 상용 근로자 데이터 Processing

Data = readmatrix("Exper_4_file.xlsx", "Range", "A2:AK145");

% 데이터 가공
dates = Data(:, 1);
sectorIndices = [2 3; 4 5; 6 7; 8 9; 10 11; 12 13; 14 15; 16 17; 18 19; 20 21; ...
                 22 23; 24 25; 26 27; 28 29; 30 31; 32 33; 34 35; 36 37];

% 폴더 생성
folderName = 'data';
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

% CSV 파일로 저장
fileNames = ["Total_data.csv", "Mine_data.csv", "Manufacturing_data.csv", "Electric_data.csv", ...
             "Waterworks_data.csv", "Construction_data.csv", "Wholesale_data.csv", ...
             "Transportation_data.csv", "Accomodation_data.csv", "IT_data.csv", ...
             "Finance_data.csv", "Real Estate_data.csv", "Scientific Service_data.csv", ...
             "Business Facility management_data.csv", "Education Service_data.csv", ...
             "Health and Social welfare services_data.csv", "Arts Sports Service_data.csv", ...
             "Repair and Other personal services_data.csv"];

for i = 1:numel(fileNames)
    sectorData = [dates, Data(:, sectorIndices(i,2)), Data(:, sectorIndices(i,1))];
    % 소수점 자리수를 조절하여 반올림 막기
    writematrix(sectorData, fileNames(i));
end
%% 일용 근로자 데이터 Processing

Data = readmatrix("Exper_5_file.xlsx", "Range", "AO4:BY147");

% 데이터 가공
dates = Data(:, 1);
sectorIndices = [2 3; 4 5; 6 7; 8 9; 10 11; 12 13; 14 15; 16 17; 18 19; 20 21; ...
                 22 23; 24 25; 26 27; 28 29; 30 31; 32 33; 34 35; 36 37];

% 폴더 생성
folderName = 'data';
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

% CSV 파일로 저장
fileNames = ["Total_data.csv", "Mine_data.csv", "Manufacturing_data.csv", "Electric_data.csv", ...
             "Waterworks_data.csv", "Construction_data.csv", "Wholesale_data.csv", ...
             "Transportation_data.csv", "Accomodation_data.csv", "IT_data.csv", ...
             "Finance_data.csv", "Real Estate_data.csv", "Scientific Service_data.csv", ...
             "Business Facility management_data.csv", "Education Service_data.csv", ...
             "Health and Social welfare services_data.csv", "Arts Sports Service_data.csv", ...
             "Repair and Other personal services_data.csv"];

for i = 1:numel(fileNames)
    sectorData = [dates, Data(:, sectorIndices(i,2)), Data(:, sectorIndices(i,1))];
    csvwrite(fileNames(i), sectorData);
end




