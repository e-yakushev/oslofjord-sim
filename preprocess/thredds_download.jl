using Dates
using Downloads

function download_thredds_files(
    start_date::Date,
    end_date::Date,
    base_url::String = "https://thredds.met.no/thredds/fileServer/fou-hi/norkyst800m-1h"
)
    current_date = start_date

    while current_date <= end_date

        # File name
        filename = "NorKyst-800m_ZDEPTHS_his.an.$(Dates.format(current_date, "yyyymmdd"))00.nc"

        # Correct URL — NO date subfolders
        url = "$base_url/$filename"

        println("Downloading $filename...")

        try
            Downloads.download(url, filename)
            println("✓ Success: $filename")

        catch e
            println("✗ Failed to download $filename: $e")
            open("failed_downloads.txt", "a") do f
                write(f, "$url\n")
            end
        end

        sleep(2)
        current_date += Day(1)
    end
end

start_date = Date(2020, 12, 16)
end_date = Date(2020, 12, 31)

download_thredds_files(start_date, end_date)

println("Done!")
